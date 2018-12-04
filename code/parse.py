import sys, os

current_line = []

for line in sys.stdin:

    split = line.split()

    if line.startswith("Input file"):

        # Print previous
        if current_line:
            print("\t".join(current_line))
            current_line = []

        current_line.append(os.path.basename(split[3][:-1]))    # Filename
        current_line.append(split[-5][:-1])                     # Number of threads
        current_line.append(split[-1])                          # Number of blocks

    elif split[-1] == "seconds":
        current_line.append(split[-2])

if current_line:
    print("\t".join(current_line))
