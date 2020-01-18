
import sys

fastfile = {}
with open("ex.txt") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            sequencenumber = line[1:]
            if sequencenumber not in fastfile:
                fastfile[sequencenumber] = []
            continue
            sequence = line
            fastfile[sequencenumber].append(sequence)

            output = []
            for i in fastfile:
                if len(fastfile[i].value()) >= 20:
                    output.append(fastfile[i].value())
                with open('output.txt', 'w') as out:
                    for i in output:
                        output.write(i)
                        output.write('\n')
                    out.close()
