
import csv
import sys

def add_unknown_lines_EDACC_csv(filename, info=False):
    all_instance_name = set()
    
    solver_inst = dict()
    with open(filename,"r") as f : 
        reader = csv.DictReader(f)
        for row in reader:
            solvername = row["Solver Configuration"]
            if solvername not in solver_inst.keys() :
                solver_inst[solvername] = list()
            solver_inst[solvername].append(row["Instance"])
            all_instance_name.add(row["Instance"])

    list_of_column_names = list()
    with open(filename,'r') as f :
        csv_reader = csv.DictReader(f)
        # converting the file to dictionary
        # by first converting to list
        # and then converting the list to dict
        dict_from_csv = dict(list(csv_reader)[0])
     
        # making a list from the keys of the dict
        list_of_column_names = list(dict_from_csv.keys())

    
    with open(filename,'a') as f :
        csv_writer = csv.DictWriter(f,fieldnames=list_of_column_names)

        for s in solver_inst.keys():
            for inst in all_instance_name :
                if inst not in solver_inst[s] :
                    dic_line = {k : "0" for k in list_of_column_names}
                    dic_line["Solver Configuration"] = s
                    dic_line["Instance"] = inst
                    dic_line["Result Code"] = "unknown"
                    dic_line["Wall Time"] = "0"
                    csv_writer.writerow(dic_line)
                    print("write line: ",inst, " for ",s,"....")
                
                
filename = sys.argv[1]
    
add_unknown_lines_EDACC_csv(filename)


    
