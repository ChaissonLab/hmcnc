import numpy as np
from scipy.special import logsumexp
from typing import Tuple, Dict, List
import warnings

class StableBaumWelch:
    def __init__(self, n_states: int):
        self.n_states = n_states
        
    def forward_pass(self, 
                    log_emissions: np.ndarray,
                    log_transitions: np.ndarray,
                    log_initial: np.ndarray) -> Tuple[np.ndarray, float]:
        """
        Compute forward probabilities in log space
        
        Args:
            log_emissions: (T, N) log emission probabilities
            log_transitions: (N, N) log transition matrix
            log_initial: (N,) log initial probabilities
            
        Returns:
            log_alpha: (T, N) log forward probabilities
            log_likelihood: float, log likelihood of sequence
        """
        T = len(log_emissions)
        log_alpha = np.zeros((T, self.n_states))
        
        # Initialize
        log_alpha[0] = log_initial + log_emissions[0]
        
        # Forward pass with scaling
        for t in range(1, T):
            # For each state j, compute log sum of incoming transitions
            for j in range(self.n_states):
                # log(sum(exp(log(alpha[t-1]) + log(trans[:,j]))))
                log_alpha[t,j] = logsumexp(
                    log_alpha[t-1] + log_transitions[:,j]
                ) + log_emissions[t,j]
                
        # Compute log likelihood using final forward values
        log_likelihood = logsumexp(log_alpha[-1])
        
        return log_alpha, log_likelihood
    
    def backward_pass(self,
                     log_emissions: np.ndarray,
                     log_transitions: np.ndarray) -> np.ndarray:
        """
        Compute backward probabilities in log space
        """
        T = len(log_emissions)
        log_beta = np.zeros((T, self.n_states))
        
        # Initialize final values to 0 (log(1) = 0)
        # Work backwards
        for t in range(T-2, -1, -1):
            for i in range(self.n_states):
                # log(sum(exp(log(trans[i,j]) + log(emit[t+1,j]) + log(beta[t+1,j]))))
                log_beta[t,i] = logsumexp(
                    log_transitions[i,:] + 
                    log_emissions[t+1,:] + 
                    log_beta[t+1,:]
                )
                
        return log_beta
    
    def compute_posteriors(self,
                         log_alpha: np.ndarray,
                         log_beta: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute state and transition posteriors
        """
        T = len(log_alpha)
        
        # State posteriors: gamma[t,i] = P(z_t = i | x)
        log_gamma = log_alpha + log_beta
        # Normalize each time step
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
        
        return gamma, xi
    
    def check_numerical_stability(self,
                                log_values: np.ndarray,
                                name: str) -> None:
        """Check for numerical issues in log-space values"""
        if np.any(np.isnan(log_values)):
            raise ValueError(f"NaN values detected in {name}")
        if np.any(np.isinf(log_values)):
            raise ValueError(f"Infinite values detected in {name}")
        
    def update_parameters(self,
                        gamma: np.ndarray,
                        xi: np.ndarray,
                        observations: np.ndarray) -> Dict:
        """
        Update model parameters using computed posteriors
        Returns dictionary of new parameter values
        """
        # Update transition probabilities
        transitions = xi.sum(axis=0)
        transitions = transitions / transitions.sum(axis=1, keepdims=True)
        
        # Update emission parameters
        # Example for negative binomial emissions:
        new_params = self._update_nb_parameters(gamma, observations)
        
        return {
            'transitions': transitions,
            'emission_params': new_params
        }
    
    def _update_nb_parameters(self,
                            gamma: np.ndarray,
                            observations: np.ndarray) -> List[Tuple[float, float]]:
        """
        Update negative binomial parameters for each state
        using weighted maximum likelihood
        """
        params = []
        for state in range(self.n_states):
            weights = gamma[:,state]
            weighted_mean = np.average(observations, weights=weights)
            weighted_var = np.average((observations - weighted_mean)**2, weights=weights)
            
            # Convert to r, p parameters
            p = weighted_mean / weighted_var
            r = weighted_mean * p / (1 - p)
            
            params.append((r, p))
            
        return params
    
    def run_estimation(self,
                      observations: np.ndarray,
                      max_iter: int = 100,
                      tol: float = 1e-6) -> Tuple[Dict, List[float]]:
        """
        Run Baum-Welch algorithm with stability checks
        """
        log_likelihoods = []
        
        # Initialize in log space
        log_transitions = np.log(self.transitions)
        log_initial = np.log(self.initial_probs)
        
        for iteration in range(max_iter):
            # Compute emission probabilities
            log_emissions = self._compute_log_emissions(observations)
            
            # Check numerical stability
            self.check_numerical_stability(log_emissions, "emissions")
            self.check_numerical_stability(log_transitions, "transitions")
            
            # Forward-backward
            try:
                log_alpha, log_likelihood = self.forward_pass(
                    log_emissions, log_transitions, log_initial
                )
                log_beta = self.backward_pass(log_emissions, log_transitions)
                
                # Store likelihood
                log_likelihoods.append(log_likelihood)
                
                # Check convergence
                if iteration > 0:
                    likelihood_change = log_likelihood - log_likelihoods[-2]
                    if abs(likelihood_change) < tol:
                        break
                        
                # Compute posteriors and update parameters
                gamma, xi = self.compute_posteriors(log_alpha, log_beta)
                new_params = self.update_parameters(gamma, xi, observations)
                
                # Convert back to log space
                log_transitions = np.log(new_params['transitions'])
                
            except (ValueError, RuntimeWarning) as e:
                warnings.warn(f"Numerical issues in iteration {iteration}: {str(e)}")
                # Could implement recovery strategies here
                break
                
        return new_params, log_likelihoods
    
    def _compute_log_emissions(self, observations: np.ndarray) -> np.ndarray:
        """
        Compute log emission probabilities for all observations and states
        """
        T = len(observations)
        log_emissions = np.zeros((T, self.n_states))
        
        for state in range(self.n_states):
            # Compute log probability for each observation under this state's distribution
            r, p = self.emission_params[state]
            log_emissions[:,state] = self._log_nb_pdf(observations, r, p)
            
        return log_emissions
    
    @staticmethod
    def _log_nb_pdf(x: np.ndarray, r: float, p: float) -> np.ndarray:
        """
        Compute log probability for negative binomial in numerically stable way
        """
        # Use logarithm of gamma functions instead of direct computation
        from scipy.special import gammaln
        
        log_prob = gammaln(x + r) - gammaln(r) - gammaln(x + 1) + \
                   r * np.log(p) + x * np.log(1 - p)
        return log_prob