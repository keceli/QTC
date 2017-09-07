
/** WORKFLOW.SWIFT
 *  Runs QTC on each line of the given input file.
 * */

import files;
import io;
import sys;
import string;

// Retrieve user argument --input=<input_list_file>
list_name = argv("input");
printf("input: " + list_name);

// Read the input file for molecule strings
list_file = input(list_name);
string lines[] = file_lines(list_file, comment="!");

// Define the app function for QTC
app qtc(string molecule)
{
  "/home/wozniak/proj/qtc/swift/qtc.sh" molecule ;
}

// Run QTC on each molecule in the input file
foreach line in lines
{
  trace(line);
  qtc(line);
}
