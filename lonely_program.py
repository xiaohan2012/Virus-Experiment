from task_dispatcher import fetch_and_do

from getopt import getopt
import sys

if __name__ == "__main__":
    opt , arg = getopt(sys.argv[1:],'n:')
    o,node = opt[0]
    fetch_and_do(node)

