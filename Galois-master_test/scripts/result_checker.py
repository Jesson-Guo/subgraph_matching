#### Script to check the output of algorithms:
### Author: Gurbinder Gill (gurbinder533@gmail.com)
### Author: Roshan Dathathri (roshan@cs.utexas.edu)
### Modified to calculate error + take tolerance as an error by Loc Hoang 

### python script.py masterFile allfile* [-t, -tolerance]=<tolerance>

### expects files in the follwing format:
###### nodeID nodeFieldVal
######## These are generated by galois::runtime::printOutput function.
### Requires python version > 2.7
# Can also take 2 field files TODO make more general

import sys
import argparse
import os

def check_results(masterFile, otherFiles, tolerance, 
  offset, errors, mrows, global_error_squared, num_nodes):

  with open(masterFile) as mfile, open(otherFiles) as ofile:
    mfile.seek(offset)

    for line2 in ofile:
      line1 = mfile.readline()
      offset = offset + len(line1)

      split_line1 = line1.split(' ')
      split_line2 = line2.split(' ')

      if (split_line1[0] == ''):
        print("ERROR: output longer than input")
        return (0, errors, mrows)

      while (long(split_line1[0]) < long(split_line2[0])):
        print "MISSING ROW: ", split_line1[0]
        mrows = mrows + 1
        line1 = mfile.readline()
        offset = offset + len(line1)
        split_line1 = line1.split(' ')


      # forces failure if missings rows exist
      #if mrows > 0:
      #  return (-1, errors, mrows)

      if (long(split_line1[0]) == long(split_line2[0])):
        # absolute value of difference in fields
        field_difference = abs(float(split_line1[1]) - float(split_line2[1]))

        global_error_squared += (field_difference ** 2)
        num_nodes += 1

        if (field_difference > tolerance):
          print "NOT MATCHED \n";
          print line1;
          print line2;
          errors = errors + 1;
        # TODO (Loc) make more general: deals with 2 fields in output (should
        # optimally deal with arbitrary # of fields
        elif (len(split_line1) == 3):
          field_difference2 = abs(float(split_line1[2]) - float(split_line2[2]))
          if (field_difference2 > tolerance):
            print "NOT MATCHED \n";
            print line1;
            print line2;
            errors = errors + 1;
      else:
        print "OFFSET MISMATCH: ", split_line1[0], split_line2[0]
        return (-1, errors, mrows, global_error_squared, num_nodes);

  return (offset, errors, mrows, global_error_squared, num_nodes);

def main(masterFile, allFiles_arr, tolerance, mean_tolerance):
  offset = 0
  errors = 0
  mrows = 0
  global_error_squared = 0
  num_nodes = 0

  for i in range(0 , len(allFiles_arr)):
    print allFiles_arr[i]
    print offset
    offset, errors, mrows, global_error_squared, num_nodes = check_results(masterFile, allFiles_arr[i], tolerance, offset, errors, mrows, global_error_squared, num_nodes)
    if (offset == -1):
      break

  rmse = (global_error_squared / num_nodes) ** 0.5
  if (rmse > mean_tolerance):
    print "\nRoot mean square error (for first field): ", rmse

  if (offset != -1):
    mfile = open(masterFile)
    mfile.seek(offset)
    old_mrows=mrows
    for line in mfile:
      mrows = mrows + 1
    if mrows > old_mrows:
      mrows = mrows - old_mrows
      print "\nNo of offsets/rows missing: ", mrows

  if (offset == -1):
    print "\nOffset not correct"

  if (errors > 0):
    print "\nNo. of mismatches: ", errors

  if (errors > 0) or (offset == -1) or (mrows > 0) or (rmse > mean_tolerance):
    print "\nFAILED\n"
    return 1
  else:
    print "\nSUCCESS\n"
    return 0

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Check graph output results")

  # parse files and an optional tolerance
  parser.add_argument('files', type=str, nargs='+', help='input + output files')
  parser.add_argument('-tolerance', '-t', type=float, nargs=1, default=0.0001,
                      help='tolerance for difference in fields (error)')
  parser.add_argument('-sort', '-s', type=bool, nargs=1, default=False,
                      help='sort the generated output files')
  parser.add_argument('-delete', '-d', type=bool, nargs=1, default=False,
                      help='delete the generated output files')
  parser.add_argument('-mean_tolerance', '-m', type=float, nargs=1, default=0.0001,
                      help='tolerance for root mean square error')

  arg = sys.argv
  parsed_arguments = parser.parse_args()

  masterFile = parsed_arguments.files[0]
  allFiles_arr = parsed_arguments.files[1:]

  print masterFile  
  print allFiles_arr  

  if parsed_arguments.sort:
    sortstr = "sort -nu"  
    for f in allFiles_arr:
      sortstr += " " + f
    sortstr += " -o .output_log"
    os.system(sortstr)

  if parsed_arguments.delete:
    rmstr = "rm -f"
    for f in allFiles_arr:
      rmstr += " " + f
    os.system(rmstr)

  if parsed_arguments.sort:
    allFiles_arr = ['.output_log']

  tolerance = parsed_arguments.tolerance
  mean_tolerance = parsed_arguments.mean_tolerance

  print("Starting comparison...")
  ret = main(masterFile, allFiles_arr, tolerance, mean_tolerance)

  if parsed_arguments.sort:
    os.system("rm -f .output_log")

  if ret:
    sys.exit(1)