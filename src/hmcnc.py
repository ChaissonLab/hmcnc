import numpy as np
from scipy.stats import nbinom
from scipy.special import logit, expit

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
    
    def emission_probability(self, coverage, clips, state):
        """
        Calculate emission probability for a window
        P(observation|state) = P(coverage|state) * P(clips|state)
        """
        # Coverage probability from negative binomial
        r, p = self.coverage_params[state]
        coverage_prob = nbinom.pmf(coverage, r, p)
        
        # Clip probability from zero-inflated negative binomial
        if clips == 0:
            clip_prob = self.clip_zero_prob + (1 - self.clip_zero_prob) * \
                       nbinom.pmf(0, self.clip_r, self.clip_p)
        else:
            clip_prob = (1 - self.clip_zero_prob) * nbinom.pmf(clips, self.clip_r, self.clip_p)
            
        return coverage_prob * clip_prob
    
    def check_convergence(self, old_params, new_params, threshold=1e-6):
        """
        Check convergence of parameter estimates
        Returns: bool indicating if converged
        """
        # Extract parameters
        old_trans, old_cov, old_clip = old_params
        new_trans, new_cov, new_clip = new_params
        
        # Calculate relative changes
        trans_change = np.max(np.abs(new_trans - old_trans) / np.abs(old_trans))
        cov_change = np.max([np.abs(new_cov[i][0] - old_cov[i][0]) / np.abs(old_cov[i][0]) 
                           for i in range(len(old_cov))])
        clip_change = np.abs(new_clip[0] - old_clip[0]) / np.abs(old_clip[0])
        
        max_change = max(trans_change, cov_change, clip_change)
        
        return max_change < threshold
    
    def baum_welch_step(self, observations):
        """
        Single step of Baum-Welch algorithm
        Returns updated parameters and log-likelihood
        """
        # Forward-backward passes would go here
        # Update equations for parameters:
        
        # 1. Update transition probabilities
        # transitions[i,j] = expected_transitions[i,j] / sum(expected_transitions[i,:])
        
        # 2. Update coverage parameters for each state
        # Use weighted MLE for negative binomial
        
        # 3. Update clipping parameters
        # Use weighted MLE for zero-inflated negative binomial
        
        # Return new parameters and log-likelihood
        pass

def run_parameter_estimation(data, max_iter=100, threshold=1e-6):
    """
    Estimate parameters using Baum-Welch with convergence checking
    """
    cnv_detector = CNVDetector()
    cnv_detector.initialize_parameters(data['coverage'], data['clips'])
    
    converged = False
    iteration = 0
    log_likelihoods = []
    
    while not converged and iteration < max_iter:
        # Store old parameters
        old_params = (
            cnv_detector.transitions.copy(),
            cnv_detector.coverage_params.copy(),
            (cnv_detector.clip_r, cnv_detector.clip_p)
        )
        
        # Run one step of Baum-Welch
        new_params, log_likelihood = cnv_detector.baum_welch_step(data)
        log_likelihoods.append(log_likelihood)
        
        # Check convergence
        converged = cnv_detector.check_convergence(old_params, new_params, threshold)
        iteration += 1
        
        # Update parameters
        cnv_detector.transitions = new_params[0]
        cnv_detector.coverage_params = new_params[1]
        cnv_detector.clip_r, cnv_detector.clip_p = new_params[2]
    
    return cnv_detector, log_likelihoods