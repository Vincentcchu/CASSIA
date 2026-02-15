#!/usr/bin/env python3
"""
Test script to verify token tracking is working correctly.
"""

import CASSIA
from CASSIA.token_tracker import get_tracker, reset_tracker

# Reset tracker
reset_tracker()
tracker = get_tracker()

# Start timer
tracker.start_timer()

# Simulate preprocessing
import time
time.sleep(0.1)
tracker.record_preprocessing_time(0.1)

# Simulate annotation
tracker.record_annotation_start()
time.sleep(0.1)
tracker.record_annotation_end()

# Simulate API calls
tracker.record_api_call(
    input_tokens=1000,
    output_tokens=500,
    model="gpt-5.1",
    provider="openai",
    purpose="annotation"
)

tracker.record_api_call(
    input_tokens=800,
    output_tokens=300,
    model="gpt-5.1",
    provider="openai",
    purpose="scoring"
)

# Stop timer
tracker.stop_timer()

# Print summary
print("\n=== Test Summary ===")
tracker.print_summary(
    input_price_per_million=1.25,
    output_price_per_million=10.0
)

# Save to file
tracker.save_summary(
    "test_performance.json",
    input_price_per_million=1.25,
    output_price_per_million=10.0
)

print("✓ Token tracking test completed successfully!")
print("✓ Performance summary saved to: test_performance.json")
