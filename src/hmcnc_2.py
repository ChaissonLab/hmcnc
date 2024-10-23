import numpy as np
from scipy.stats import nbinom
from scipy.special import logit, expit, logsumexp
from typing import Dict, Tuple, List

class CNVDetector:
    def __init__(self, n_states=5):
        """
        Initialize CNV detector with states representing copy numbers 0-4
        n_states: number of copy number states (default 5 for CN 0-4)
        """
        self.n_states = n_states
        # Initial uniform state probabilities
        self.initial_probs = np.ones(n_states) / n_states
        
    def initialize_parameters(self, coverage_data, clip_data):
        """
        Initialize model parameters using data statistics
        """
        # Zero-inflated negative binomial parameters for clipping
        self.clip_zero_prob = np.sum(clip_data == 0) / len(clip_data)
        self.clip_r, self.clip_p = self._fit_nb_initial(clip_data[clip_data > 0])
        
        # Negative binomial parameters for coverage per copy number state
        coverage_mean = np.mean(coverage_data)
        coverage_var = np.var(coverage_data)
        
        # Initialize parameters for each copy number state
        self.coverage_params = []
        for cn in range(self.n_states):
            # Expected coverage scales with copy number
            expected_mean = coverage_mean * (cn / 2)  # Assuming CN=2 is normal
            # Maintain similar dispersion
            r, p = self._convert_mean_var_to_nb(expected_mean, coverage_var)
            self.coverage_params.append((r, p))
        
        # Initialize transition matrix with strong diagonal
        self.transitions = np.eye(self.n_states) * 0.9
        self.transitions += np.ones((self.n_states, self.n_states)) * 0.1 / self.n_states
        
    def _convert_mean_var_to_nb(self, mean, var):
        """Convert mean and variance to negative binomial parameters"""
        p = mean / var
        r = mean * p / (1 - p)
        return r, p
        
    def _fit_nb_initial(self, data):
        """Fit initial negative binomial parameters"""
        mean = np.mean(data)
        var = np.var(data)
        return self._convert_mean_var_to_nb(mean, var)
    
    def _compute_log_emission_probs(self, observations: Dict) -> np.ndarray:
        """
        Compute log emission probabilities for all observations and states
        Returns: (T, n_states) array of log probabilities
        """
        T = len(observations['coverage'])
        log_emissions = np.zeros((T, self.n_states))
        
        for t in range(T):
            coverage = observations['coverage'][t]
            clips = observations['clips'][t]
            
            for state in range(self.n_states):
                # Coverage probability
                r, p = self.coverage_params[state]
                log_coverage_prob = nbinom.logpmf(coverage, r, p)
                
                # Clip probability
                if clips == 0:
                    clip_prob = self.clip_zero_prob + (1 - self.clip_zero_prob) * \
                               nbinom.pmf(0, self.clip_r, self.clip_p)
                    log_clip_prob = np.log(clip_prob)
                else:
                    log_clip_prob = np.log(1 - self.clip_zero_prob) + \
                                  nbinom.logpmf(clips, self.clip_r, self.clip_p)
                
                log_emissions[t, state] = log_coverage_prob + log_clip_prob
                
        return log_emissions
    
    def baum_welch_step(self, observations: Dict) -> Tuple[tuple, float]:
        """
        Single step of Baum-Welch algorithm
        Returns:
            - Tuple of (transitions, coverage_params, clip_params)
            - Log likelihood
        """
        # Get log emission probabilities
        log_emissions = self._compute_log_emission_probs(observations)
        log_transitions = np.log(self.transitions)
        log_initial = np.log(self.initial_probs)
        
        T = len(observations['coverage'])
        
        # Forward pass
        log_alpha = np.zeros((T, self.n_states))
        log_alpha[0] = log_initial + log_emissions[0]
        
        for t in range(1, T):
            for j in range(self.n_states):
                log_alpha[t,j] = logsumexp(
                    log_alpha[t-1] + log_transitions[:,j]
                ) + log_emissions[t,j]
        
        # Calculate log likelihood
        log_likelihood = logsumexp(log_alpha[-1])
        
        # Backward pass
        log_beta = np.zeros((T, self.n_states))
        for t in range(T-2, -1, -1):
            for i in range(self.n_states):
                log_beta[t,i] = logsumexp(
                    log_transitions[i,:] + 
                    log_emissions[t+1,:] + 
                    log_beta[t+1,:]
                )
        
        # Compute posteriors
        # State posteriors: gamma[t,i] = P(z_t = i | x)
        log_gamma = log_alpha + log_beta
        # Normalize
        for t in range(T):
            log_gamma[t] -= logsumexp(log_gamma[t])
        gamma = np.exp(log_gamma)
        
        # Transition posteriors: xi[t,i,j] = P(z_t = i, z_{t+1} = j | x)
        xi = np.zeros((T-1, self.n_states, self.n_states))
        for t in range(T-1):
            for i in range(self.n_states):
                for j in range(self.n_states):
                    xi[t,i,j] = np.exp(
                        log_alpha[t,i] +
                        log_transitions[i,j] +
                        log_emissions[t+1,j] +
                        log_beta[t+1,j] -
                        log_likelihood
                    )
        
        # Update parameters
        # 1. Transition probabilities
        new_transitions = xi.sum(axis=0)
        new_transitions /= new_transitions.sum(axis=1, keepdims=True)
        
        # 2. Coverage parameters for each state
        new_coverage_params = []
        for state in range(self.n_states):
            weights = gamma[:, state]
            weighted_mean = np.average(observations['coverage'], weights=weights)
            weighted_var = np.average(
                (observations['coverage'] - weighted_mean)**2, 
                weights=weights
            )
            r, p = self._convert_mean_var_to_nb(weighted_mean, weighted_var)
            new_coverage_params.append((r, p))
        
        # 3. Clipping parameters
        # Update zero probability
        new_clip_zero_prob = np.average(
            observations['clips'] == 0,
            weights=gamma.sum(axis=1)
        )
        
        # Update negative binomial parameters for non-zero counts
        non_zero_mask = observations['clips'] > 0
        if np.any(non_zero_mask):
            weights = gamma.sum(axis=1)[non_zero_mask]
            weighted_mean = np.average(
                observations['clips'][non_zero_mask], 
                weights=weights
            )
            weighted_var = np.average(
                (observations['clips'][non_zero_mask] - weighted_mean)**2,
                weights=weights
            )
            new_clip_r, new_clip_p = self._convert_mean_var_to_nb(
                weighted_mean, weighted_var
            )
        else:
            new_clip_r, new_clip_p = self.clip_r, self.clip_p
        
        return (
            (new_transitions, new_coverage_params, (new_clip_r, new_clip_p)),
            log_likelihood
        )

    def viterbi(self, observations: Dict) -> List[int]:
        """
        Run Viterbi algorithm to find most likely state sequence
        """
        T = len(observations['coverage'])
        log_emissions = self._compute_log_emission_probs(observations)
        log_transitions = np.log(self.transitions)
        log_initial = np.log(self.initial_probs)
        
        # Initialize
        viterbi_matrix = np.zeros((T, self.n_states))
        backpointer = np.zeros((T, self.n_states), dtype=int)
        
        # Initial step
        viterbi_matrix[0] = log_initial + log_emissions[0]
        
        # Recursion
        for t in range(1, T):
            for j in range(self.n_states):
                prob_vector = viterbi_matrix[t-1] + log_transitions[:,j]
                backpointer[t,j] = np.argmax(prob_vector)
                viterbi_matrix[t,j] = np.max(prob_vector) + log_emissions[t,j]
        
        # Backtrack
        states = np.zeros(T, dtype=int)
        states[-1] = np.argmax(viterbi_matrix[-1])
        for t in range(T-2, -1, -1):
            states[t] = backpointer[t+1, states[t+1]]
        
        return states.tolist()