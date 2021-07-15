import os
import sys
import subprocess

# Remesh python script
# will take filename or relative path to look into
def main():
    filename = sys.argv[1]
    log = open(filename)
    last_lines = tail (log,lines=100)
    log.close()
    for line in last_lines:
        if "Saturated angle" in line:
            angle = line.split()[2].replace(';','')
            print("Found deflection angle= ", angle)
    
    replace_line("parameters.dat","rotateAirfoil="+angle+'\n',2)
    print("Calling mesher")
    subprocess.run(sys.argv[2:])

def usage():
    print("----------SCRIPT USAGE------------")
    print("python3 logfile mesherName mesherOptions")
    print("----------------------------------")


def replace_line(filename,new,pos):
    with open(filename,'r') as f:
        filelines = f.readlines()
    
    filelines[pos] = new

    with open(filename,'w') as f:
        f.writelines(filelines)
    
def tail(f, lines=1, _buffer=4098):
    """Tail a file and get X lines from the end"""
    # place holder for the lines found
    lines_found = []

    # block counter will be multiplied by buffer
    # to get the block size from the end
    block_counter = -1

    # loop until we find X lines
    while len(lines_found) < lines:
        try:
            f.seek(block_counter * _buffer, os.SEEK_END)
        except IOError:  # either file is too small, or too many lines requested
            f.seek(0)
            lines_found = f.readlines()
            break

        lines_found = f.readlines()

        # we found enough lines, get out
        # Removed this line because it was redundant the while will catch
        # it, I left it for history
        # if len(lines_found) > lines:
        #    break

        # decrement the block counter to get the
        # next X bytes
        block_counter -= 1

    return lines_found[-lines:]

if __name__ == "__main__":
    usage()
    main()
else:
    print("DO NOTHING. SCRIPT ONLY CALLABLE FROM CMDLINE")