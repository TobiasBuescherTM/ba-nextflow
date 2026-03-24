import os
from datetime import datetime
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--log_file", required=True )
parser.add_argument("--log_output", required=True )


args = parser.parse_args()

now = datetime.now()
timestamp = now.strftime("%Y-%m-%d %H:%M:%S")

with open(args.log_file, "r") as src, open(args.log_output, "w") as dst:
    dst.write(src.read())

with open(args.log_output, "a") as file:
    file.write(f"vsnStart: {timestamp}\n")