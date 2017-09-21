
/** WORKFLOW.SWIFT
 *  Runs QTC on each line of the given input file.
 * */

import files;
import io;
import string;
import sys;

// Project root directory
string THIS = getenv("THIS"); 

// Retrieve user argument --input=<input_list_file>
list_name = argv("input");
printf("input: " + list_name);

// Read the input file for molecule strings
list_file = input(list_name);
string lines[] = file_lines(list_file, comment="!");

// Define the app function for QTC
app qtc(string molecule, int index)
{
  (THIS+"/qtc.sh") molecule index ;
}

// Run QTC on each molecule in the input file
foreach line,i in lines
{
  trace(line + " " + i);
  // trace(hash(line));
  qtc(line, i);
}
