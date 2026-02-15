"""
Token and cost tracking module for CASSIA.

This module tracks API calls, token usage, and costs across the CASSIA pipeline.
"""

import time
import json
from typing import Dict, Any, Optional
from threading import RLock

class TokenTracker:
    """Thread-safe token and cost tracker for LLM API calls."""
    
    def __init__(self):
        self.lock = RLock()  # Use RLock (reentrant lock) to allow nested lock acquisition
        self.reset()
    
    def reset(self):
        """Reset all tracking counters."""
        with self.lock:
            self.total_input_tokens = 0
            self.total_output_tokens = 0
            self.total_calls = 0
            self.start_time = None
            self.preprocessing_time = 0
            self.annotation_time = 0
            self.total_time = 0
            self.call_details = []  # List of individual call details
    
    def start_timer(self):
        """Start the overall timer."""
        with self.lock:
            self.start_time = time.time()
    
    def record_preprocessing_time(self, duration: float):
        """Record preprocessing time in seconds."""
        with self.lock:
            self.preprocessing_time = duration
    
    def record_annotation_start(self):
        """Mark the start of annotation phase."""
        with self.lock:
            self.annotation_start_time = time.time()
    
    def record_annotation_end(self):
        """Mark the end of annotation phase."""
        with self.lock:
            if hasattr(self, 'annotation_start_time'):
                self.annotation_time = time.time() - self.annotation_start_time
    
    def stop_timer(self):
        """Stop the overall timer."""
        with self.lock:
            if self.start_time:
                self.total_time = time.time() - self.start_time
    
    def record_api_call(
        self,
        input_tokens: int,
        output_tokens: int,
        model: str,
        provider: str,
        purpose: str = "annotation"
    ):
        """
        Record an API call with token usage.
        
        Args:
            input_tokens: Number of input tokens used
            output_tokens: Number of output tokens used
            model: Model name used
            provider: Provider name
            purpose: Purpose of the call (e.g., "annotation", "scoring", "boost")
        """
        with self.lock:
            self.total_input_tokens += input_tokens
            self.total_output_tokens += output_tokens
            self.total_calls += 1
            self.call_details.append({
                'input_tokens': input_tokens,
                'output_tokens': output_tokens,
                'model': model,
                'provider': provider,
                'purpose': purpose
            })
    
    def calculate_cost(
        self,
        input_price_per_million: float = 1.25,
        output_price_per_million: float = 10.0
    ) -> float:
        """
        Calculate total cost based on token usage.
        
        Args:
            input_price_per_million: Price per 1M input tokens (default: $1.25 for GPT-5.1)
            output_price_per_million: Price per 1M output tokens (default: $10 for GPT-5.1)
        
        Returns:
            Total cost in dollars
        """
        with self.lock:
            input_cost = (self.total_input_tokens / 1_000_000) * input_price_per_million
            output_cost = (self.total_output_tokens / 1_000_000) * output_price_per_million
            return input_cost + output_cost
    
    def get_summary(
        self,
        input_price_per_million: float = 1.25,
        output_price_per_million: float = 10.0
    ) -> Dict[str, Any]:
        """
        Get a summary of all tracking data.
        
        Args:
            input_price_per_million: Price per 1M input tokens
            output_price_per_million: Price per 1M output tokens
        
        Returns:
            Dictionary with summary statistics
        """
        with self.lock:
            total_cost = self.calculate_cost(input_price_per_million, output_price_per_million)
            
            return {
                'timing': {
                    'preprocessing_seconds': round(self.preprocessing_time, 2),
                    'annotation_seconds': round(self.annotation_time, 2),
                    'total_seconds': round(self.total_time, 2),
                    'preprocessing_minutes': round(self.preprocessing_time / 60, 2),
                    'annotation_minutes': round(self.annotation_time / 60, 2),
                    'total_minutes': round(self.total_time / 60, 2)
                },
                'tokens': {
                    'total_input_tokens': self.total_input_tokens,
                    'total_output_tokens': self.total_output_tokens,
                    'total_tokens': self.total_input_tokens + self.total_output_tokens
                },
                'api_calls': {
                    'total_calls': self.total_calls
                },
                'cost': {
                    'input_cost_usd': round((self.total_input_tokens / 1_000_000) * input_price_per_million, 4),
                    'output_cost_usd': round((self.total_output_tokens / 1_000_000) * output_price_per_million, 4),
                    'total_cost_usd': round(total_cost, 4),
                    'input_price_per_million': input_price_per_million,
                    'output_price_per_million': output_price_per_million
                }
            }
    
    def print_summary(
        self,
        input_price_per_million: float = 1.25,
        output_price_per_million: float = 10.0
    ):
        """Print a formatted summary of tracking data."""
        summary = self.get_summary(input_price_per_million, output_price_per_million)
        
        print("\n" + "="*70)
        print("CASSIA Pipeline Performance Summary")
        print("="*70)
        
        print("\n⏱️  Timing:")
        print(f"  Preprocessing: {summary['timing']['preprocessing_minutes']:.2f} min")
        print(f"  Annotation:    {summary['timing']['annotation_minutes']:.2f} min")
        print(f"  Total:         {summary['timing']['total_minutes']:.2f} min")
        
        print("\n🎯 Token Usage:")
        print(f"  Input tokens:  {summary['tokens']['total_input_tokens']:,}")
        print(f"  Output tokens: {summary['tokens']['total_output_tokens']:,}")
        print(f"  Total tokens:  {summary['tokens']['total_tokens']:,}")
        
        print("\n📞 API Calls:")
        print(f"  Total calls:   {summary['api_calls']['total_calls']}")
        
        print("\n💰 Cost (USD):")
        print(f"  Input cost:    ${summary['cost']['input_cost_usd']:.4f}")
        print(f"  Output cost:   ${summary['cost']['output_cost_usd']:.4f}")
        print(f"  Total cost:    ${summary['cost']['total_cost_usd']:.4f}")
        print(f"  (Rates: ${input_price_per_million}/M input, ${output_price_per_million}/M output)")
        
        print("="*70 + "\n")
    
    def save_summary(
        self,
        filepath: str,
        input_price_per_million: float = 1.25,
        output_price_per_million: float = 10.0
    ):
        """Save summary to JSON file."""
        summary = self.get_summary(input_price_per_million, output_price_per_million)
        with open(filepath, 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"Performance summary saved to: {filepath}")


# Global tracker instance
_global_tracker = TokenTracker()

def get_tracker() -> TokenTracker:
    """Get the global token tracker instance."""
    return _global_tracker

def reset_tracker():
    """Reset the global tracker."""
    _global_tracker.reset()
