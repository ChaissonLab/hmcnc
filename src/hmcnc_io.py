import pandas as pd
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import List, Tuple, Dict
import logging

@dataclass
class WindowData:
    """Container for 100bp window data"""
    chrom: str
    start: int
    end: int
    coverage: int
    left_clips: int
    right_clips: int
    
    @property
    def total_clips(self) -> int:
        return self.left_clips + self.right_clips

class CNVPipeline:
    def __init__(self, window_size: int = 100):
        self.window_size = window_size
        self.windows: Dict[str, List[WindowData]] = {}
        self.detector = None
        
    def load_bed_files(self, 
                      coverage_bed: Path,
                      left_clips_bed: Path,
                      right_clips_bed: Path) -> None:
        """
        Load and merge coverage and clipping BED files
        
        Format expected:
        chrom start end count
        """
        # Load files
        coverage_df = pd.read_csv(coverage_bed, sep='\t', 
                                names=['chrom', 'start', 'end', 'coverage'])
        left_df = pd.read_csv(left_clips_bed, sep='\t',
                            names=['chrom', 'start', 'end', 'left_clips'])
        right_df = pd.read_csv(right_clips_bed, sep='\t',
                             names=['chrom', 'start', 'end', 'right_clips'])
        
        # Merge on genomic coordinates
        merged = coverage_df.merge(left_df, on=['chrom', 'start', 'end'])
        merged = merged.merge(right_df, on=['chrom', 'start', 'end'])
        
        # Convert to WindowData objects
        self.windows = {
            chrom: [
                WindowData(row.chrom, row.start, row.end, 
                         row.coverage, row.left_clips, row.right_clips)
                for _, row in group.iterrows()
            ]
            for chrom, group in merged.groupby('chrom')
        }
        
        logging.info(f"Loaded {sum(len(wins) for wins in self.windows.values())} windows")
        
    def prepare_training_data(self, 
                            chrom: str,
                            exclude_regions: List[Tuple[int, int]] = None) -> Dict:
        """
        Prepare training data for parameter estimation
        
        Args:
            chrom: Chromosome to use for training
            exclude_regions: List of (start, end) tuples to exclude
        """
        windows = self.windows[chrom]
        
        if exclude_regions:
            windows = [
                w for w in windows
                if not any(start <= w.start <= end 
                          for start, end in exclude_regions)
            ]
        
        return {
            'coverage': np.array([w.coverage for w in windows]),
            'clips': np.array([w.total_clips for w in windows])
        }
    
    def segment_chromosome(self, 
                         chrom: str,
                         min_segment_size: int = 5) -> List[Dict]:
        """
        Call CNV segments for a chromosome
        
        Returns list of segments:
        [{
            'chrom': str,
            'start': int,
            'end': int,
            'copy_number': int,
            'confidence': float
        }]
        """
        windows = self.windows[chrom]
        
        # Prepare observation sequence
        observations = {
            'coverage': np.array([w.coverage for w in windows]),
            'clips': np.array([w.total_clips for w in windows])
        }
        
        # Run Viterbi algorithm to get most likely state sequence
        states = self.detector.viterbi(observations)
        
        # Merge adjacent windows with same state
        segments = []
        current_segment = {
            'chrom': chrom,
            'start': windows[0].start,
            'state': states[0],
            'windows': [windows[0]]
        }
        
        for window, state in zip(windows[1:], states[1:]):
            if state == current_segment['state']:
                current_segment['windows'].append(window)
            else:
                # Finish current segment
                if len(current_segment['windows']) >= min_segment_size:
                    segments.append({
                        'chrom': chrom,
                        'start': current_segment['start'],
                        'end': current_segment['windows'][-1].end,
                        'copy_number': current_segment['state'],
                        'confidence': self._calculate_segment_confidence(
                            current_segment['windows'],
                            current_segment['state']
                        )
                    })
                
                # Start new segment
                current_segment = {
                    'chrom': chrom,
                    'start': window.start,
                    'state': state,
                    'windows': [window]
                }
        
        # Add final segment
        if len(current_segment['windows']) >= min_segment_size:
            segments.append({
                'chrom': chrom,
                'start': current_segment['start'],
                'end': current_segment['windows'][-1].end,
                'copy_number': current_segment['state'],
                'confidence': self._calculate_segment_confidence(
                    current_segment['windows'],
                    current_segment['state']
                )
            })
        
        return segments
    
    def _calculate_segment_confidence(self, 
                                   windows: List[WindowData], 
                                   state: int) -> float:
        """Calculate confidence score for a segment"""
        # Example confidence calculation:
        # - How well coverage matches expected for copy number
        # - Consistency of coverage across segment
        # - Support from soft clipping
        coverages = [w.coverage for w in windows]
        clips = [w.total_clips for w in windows]
        
        expected_coverage = np.mean(coverages) * (state / 2)
        coverage_deviation = np.std(coverages) / expected_coverage
        
        # Return confidence score 0-1
        return 1 / (1 + coverage_deviation)

def run_cnv_analysis(coverage_bed: Path,
                    left_clips_bed: Path,
                    right_clips_bed: Path,
                    training_chrom: str = 'chr1',
                    exclude_regions: List[Tuple[int, int]] = None) -> Dict:
    """
    Run complete CNV analysis pipeline
    
    Returns dictionary of CNV segments by chromosome
    """
    # Initialize pipeline
    pipeline = CNVPipeline()
    
    # Load data
    pipeline.load_bed_files(coverage_bed, left_clips_bed, right_clips_bed)
    
    # Prepare training data
    training_data = pipeline.prepare_training_data(
        training_chrom, exclude_regions
    )
    
    # Initialize and train CNV detector
    pipeline.detector = CNVDetector()
    pipeline.detector, log_likelihoods = run_parameter_estimation(training_data)
    
    # Call segments for each chromosome
    results = {}
    for chrom in pipeline.windows.keys():
        results[chrom] = pipeline.segment_chromosome(chrom)
    
    return results