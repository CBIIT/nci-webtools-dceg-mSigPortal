
import subprocess
import sys  
import json
if __name__ == '__main__':
    script = "test.r"
    commands = ["rscript", "api\\r\\" + script]
    arguments = sys.argv[1:]
    for args in arguments:
        json.dumps(args)

    subprocess.call(commands + arguments, shell=True)
