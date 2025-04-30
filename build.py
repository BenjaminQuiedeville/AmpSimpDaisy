import os 
import sys 

def run(command):
    print(command, file = sys.stdout, flush = True)
    result = os.system(command)
    return result

run("make clean")

result = 0

rewrite_ir = False

if rewrite_ir:
    result = run("pushd data && julia --startup-file=no ir_generation.jl && popd")
    if result != 0: exit(result)


result = run("make")

if result != 0: exit(result)

result = run("make program")

exit(result)
