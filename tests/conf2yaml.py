import os
import re

def string_clean(str_in):
    str_in=str_in.replace("  "," ").strip()
    str_in=str_in.replace("  "," ").strip()
    str_in=str_in.replace("  "," ").strip()
    str_in=str_in.replace("  "," ").strip()
    str_out="'"+str_in.replace(" ","','")+"'"
    return str_out

if __name__ == "__main__":
    with open('rt.yaml', 'w') as yaml_file, open("rt.conf") as conf_file:
        for line in conf_file:
            line = line.strip()
            if not line:  # skip: line is blank
                continue
            if line.startswith("#"):  # skip: comment line
                continue
            if line.startswith("COMPILE"):  # COMPILE line
                build    = line.split('|')
                apps     =     build[1].strip()
                compiler = "'"+build[2].strip()+"'"
                options  = "'"+build[3].strip()+"'"
                machine  =     build[4].strip()
                off_machine = None
                if (machine.find('-') != -1):
                    off_machine=machine.replace("-","").strip()
                    off_machine=string_clean(off_machine)
                yaml_file.write(apps+":"+ '\n')
                yaml_file.write("  build: "+ '\n')
                yaml_file.write("    compiler: "+compiler+ '\n')
                yaml_file.write("    option: "+options+ '\n')
                if not (off_machine is None):
                    yaml_file.write("    turnoff: ["+off_machine+"]"+ '\n')
                prev_line='COMPILE'
            if line.startswith("RUN"):  # RUN line
                build    = line.split('|')
                test     =     build[1].strip()
                machine  =     build[2].strip()
                baseline = "'"+build[3].strip()+"'"
                depend   =     build[4].strip()
                if (machine.find('-') != -1):
                    off_machine=machine.replace("-","").strip()
                    off_machine=string_clean(off_machine)
                tests = "    "+"- "+test+": {'recurring':['daily']"
                if baseline.isalnum(): tests = tests + ",'baseline': "+baseline
                if depend and depend.strip(): tests = tests + ",'dependency':'"+depend+"'"                
                if not (off_machine is None): tests = tests +",'turnoff':["+off_machine+"]"
                if prev_line == "COMPILE": yaml_file.write("  tests: "+ '\n')
                yaml_file.write(tests+"}"+ '\n')                
                prev_line='RUN'
        
    yaml_file.close(); conf_file.close()

