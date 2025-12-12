import sys
import os

print("Python executable:", sys.executable)
print("\nSearching for logging.py files...\n")

for path in sys.path:
    if os.path.isdir(path):
        for root, dirs, files in os.walk(path):
            if 'logging.py' in files:
                full_path = os.path.join(root, 'logging.py')
                print(f"FOUND: {full_path}")