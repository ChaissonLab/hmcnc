import numpy as np
from scipy.stats import nbinom
from scipy.special import logit, expit, logsumexp
from typing import Dict, Tuple, List
from enum import Enum

class CNVState(Enum):
    """Enum for CNV states with associated copy numbers"""
    DELETION = 0
    DELETION_START = 1  # Higher left clips expected
    DELETION_END = 2    # Higher right clips expected
    NORMAL = 3
    DUPLICATION_START = 4  # Higher right clips expected
    DUPLICATION_END = 5    # Higher left clips expected
    DUPLICATION = 6

class CNVDetector:
    def __init__(self, n_copy_states=4):
        """
        Initialize CNV detector with states for:
        - Each copy number (0-3 by default)
        - Boundary states (deletion start/end, duplication start/end)
        """
        self.n_copy_states = n_copy_states
        self.n_states = len(CNVState)
        # Initial state probabilities - higher for normal state
        self.initial_probs = np.ones(self.n_states) * 0.1
        self.initial_probs[CNVState.NORMAL.value] = 0.4
        self.initial_probs /= np.sum(self.initial_probs)
        
    def initialize_parameters(self, coverage_data, left_clip_data, right_clip_data):
        """
        Initialize model parameters using data statistics
        Now handles left and right clips separately
        """
        # Initialize clip parameters for each direction
        self.left_clip_zero_prob = np.sum(left_clip_data == 0) / len(left_clip_data)
        self.right_clip_zero_prob = np.sum(right_clip_data == 0) / len(right_clip_data)
        
        # Fit clip parameters for non-zero counts
        self.left_clip_r, self.left_clip_p = self._fit_nb_initial(left_clip_data[left_clip_data > 0])
        self.right_clip_r, self.right_clip_p = self._fit_nb_initial(right_clip_data[right_clip_data > 0])
        
        # Coverage parameters for each copy number
        coverage_mean = np.mean(coverage_data)
        coverage_var = np.var(coverage_data)
        
        # Initialize coverage parameters for each state
        self.coverage_params = []
        for state in CNVState:
            if state == CNVState.NORMAL:
                cn = 2
            elif state in [CNVState.DELETION, CNVState.DELETION_START, CNVState.DELETION_END]:
                cn = 1
            elif state in [CNVState.DUPLICATION, CNVState.DUPLICATION_START, CNVState.DUPLICATION_END]:
                cn = 3
            else:
                cn = 2  # Default for boundary states
                
            expected_mean = coverage_mean * (cn / 2)
            r, p = self._convert_mean_var_to_nb(expected_mean, coverage_var)
            self.coverage_params.append((r, p))
        
        # Initialize transition matrix
        self.transitions = self._initialize_transitions()
        
    def _initialize_transitions(self) -> np.ndarray:
        """
        Initialize transition matrix with biological constraints:
        - High probability of staying in same state
        - Transitions through boundary states when changing copy number
        - No direct transitions between deletion and duplication
        """
        n = self.n_states
        transitions = np.zeros((n, n))
        
        # High probability of staying in same state
        for i in range(n):
            transitions[i,i] = 0.8
        
        # Define allowed transitions
        allowed_transitions = {
            CNVState.NORMAL: [CNVState.DELETION_START, CNVState.DUPLICATION_START],
            CNVState.DELETION_START: [CNVState.DELETION],
            CNVState.DELETION: [CNVState.DELETION_END],
            CNVState.DELETION_END: [CNVState.NORMAL],
            CNVState.DUPLICATION_START: [CNVState.DUPLICATION],
            CNVState.DUPLICATION: [CNVState.DUPLICATION_END],
            CNVState.DUPLICATION_END: [CNVState.NORMAL]
        }
        
        # Set transition probabilities
        for state in CNVState:
            if state in allowed_transitions:
                prob = (1 - transitions[state.value, state.value]) / len(allowed_transitions[state])
                for next_state in allowed_transitions[state]:
                    transitions[state.value, next_state.value] = prob
        
        return transitions
    
    def _compute_log_emission_probs(self, observations: Dict) -> np.ndarray:
        """
        Compute log emission probabilities accounting for directional clipping
        """
        T = len(observations['coverage'])
        log_emissions = np.zeros((T, self.n_states))
        
        for t in range(T):
            coverage = observations['coverage'][t]
            left_clips = observations['left_clips'][t]
            right_clips = observations['right_clips'][t]
            
            for state in CNVState:
                # 1. Coverage probability
                r, p = self.coverage_params[state.value]
                log_coverage_prob = nbinom.logpmf(coverage, r, p)
                
                # 2. Clipping probabilities with directional expectations
                log_clip_prob = self._compute_directional_clip_prob(
                    state, left_clips, right_clips
                )
                
                log_emissions[t, state.value] = log_coverage_prob + log_clip_prob
                
        return log_emissions
    
    def _compute_directional_clip_prob(self, 
                                     state: CNVState, 
                                     left_clips: int, 
                                     right_clips: int) -> float:
        """
        Compute clipping probability based on state and expected clip direction
        """
        # Get base probabilities
        left_prob = self._compute_single_clip_prob(
            left_clips, self.left_clip_zero_prob,
            self.left_clip_r, self.left_clip_p
        )
        right_prob = self._compute_single_clip_prob(
            right_clips, self.right_clip_zero_prob,
            self.right_clip_r, self.right_clip_p
        )
        
        # Adjust probabilities based on state expectations
        if state == CNVState.DELETION_START:
            # Expect more left clips
            left_prob *= 2
        elif state == CNVState.DELETION_END:
            # Expect more right clips
            right_prob *= 2
        elif state == CNVState.DUPLICATION_START:
            # Expect more right clips
            right_prob *= 2
        elif state == CNVState.DUPLICATION_END:
            # Expect more left clips
            left_prob *= 2
            
        return np.log(left_prob) + np.log(right_prob)
    
    def _compute_single_clip_prob(self, 
                                clips: int, 
                                zero_prob: float,
                                r: float, 
                                p: float) -> float:
        """Compute probability for single direction of clipping"""
        if clips == 0:
            return zero_prob + (1 - zero_prob) * nbinom.pmf(0, r, p)
        else:
            return (1 - zero_prob) * nbinom.pmf(clips, r, p)
    
    def baum_welch_step(self, observations: Dict) -> Tuple[tuple, float]:
        """
        Single step of Baum-Welch algorithm with directional clip handling
        """
        # Forward-backward passes remain the same
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
        log_gamma = log_alpha + log_beta
        for t in range(T):
            log_gamma[t] -= logsumexp(log_gamma[t])
        gamma = np.exp(log_gamma)
        
        # Update parameters
        new_transitions = self._update_transitions(
            gamma, log_alpha, log_beta, log_emissions, log_transitions
        )
        
        new_coverage_params = self._update_coverage_params(
            gamma, observations['coverage']
        )
        
        new_clip_params = self._update_clip_params(
            gamma, 
            observations['left_clips'],
            observations['right_clips']
        )
        
        return (new_transitions, new_coverage_params, new_clip_params), log_likelihood
    
    def _update_clip_params(self, 
                          gamma: np.ndarray,
                          left_clips: np.ndarray,
                          right_clips: np.ndarray) -> tuple:
        """Update clipping parameters for both directions"""
        # Update zero probabilities
        new_left_zero_prob = np.average(
            left_clips == 0,
            weights=gamma.sum(axis=1)
        )
        new_right_zero_prob = np.average(
            right_clips == 0,
            weights=gamma.sum(axis=1)
        )
        
        # Update negative binomial parameters
        left_params = self._update_single_clip_params(gamma, left_clips)
        right_params = self._update_single_clip_params(gamma, right_clips)
        
        return (
            (new_left_zero_prob, *left_params),
            (new_right_zero_prob, *right_params)
        )
    
    def _update_single_clip_params(self,
                                 gamma: np.ndarray,
                                 clips: np.ndarray) -> Tuple[float, float]:
        """Update parameters for single direction of clipping"""
        non_zero_mask = clips > 0
        if np.any(non_zero_mask):
            weights = gamma.sum(axis=1)[non_zero_mask]
            weighted_mean = np.average(clips[non_zero_mask], weights=weights)
            weighted_var = np.average(
                (clips[non_zero_mask] - weighted_mean)**2,
                weights=weights
            )
            return self._convert_mean_var_to_nb(weighted_mean, weighted_var)
        return self.left_clip_r, self.left_clip_p  # Return old parameters if no non-zero data
    
    def viterbi(self, observations: Dict) -> List[CNVState]:
        """
        Run Viterbi algorithm to find most likely state sequence
        Returns list of CNVState enum values
        """
        T = len(observations['coverage'])
        log_emissions = self._compute_log_emission_probs(observations)
        log_transitions = np.log(self.transitions)
        log_initial = np.log(self.initial_probs)
        
        viterbi_matrix = np.zeros((T, self.n_states))
        backpointer = np.zeros((T, self.n_states), dtype=int)
        
        viterbi_matrix[0] = log_initial + log_emissions[0]
        
        for t in range(1, T):
            for j in range(self.n_states):
                prob_vector = viterbi_matrix[t-1] + log_transitions[:,j]
                backpointer[t,j] = np.argmax(prob_vector)
                viterbi_matrix[t,j] = np.max(prob_vector) + log_emissions[t,j]
        
        states = np.zeros(T, dtype=int)
        states[-1] = np.argmax(viterbi_matrix[-1])
        for t in range(T-2, -1, -1):
            states[t] = backpointer[t+1, states[t+1]]
        
        return [CNVState(s) for s in states]