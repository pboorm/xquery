# xquery
X-ray Query - find X-ray observations of positions on the sky

Author: Peter Boorman, Czech Academy of Sciences

## Usage:
xquery.py --textfile --wgetlocs --filterdata --columns

 --textfile [TEXTFILE], -t [TEXTFILE]
                        Name of textfile containing names and/or coordinates,
                        one per line. Default textfile name == 'mysources.txt'
                        
 --wgetlocs [WGETLOCS], -w [WGETLOCS]
                        Choice to find the wget locations of all data found.
                        Warning: takes a long time. Can be y|Y|n|N, default ==
                        n.
                        
 --filterdata [FILTERDATA], -f [FILTERDATA]
                        Choice to filter the observations found for data that
                        is archived, i.e. available. Can be y|Y|n|N, default
                        == n.
                        
 --columns [COLUMNS], -c [COLUMNS]
                        Choice to include 'standard' columns or 'all' columns
                        in returned tables.

`xquery` searches the High Energy Astrophysics Science Archive Research Center (HEASARC) archive for any archival X-ray data available for a user-provided list of targets. The publicly available Javascript (users/jar) is used in `xquery`, which is also available from here: https://heasarc.gsfc.nasa.gov/xamin/doc/CLIUsersGuide.html

- Input:
Targets can be provided in a text file, with one source identifier (or set of coordinates) per line. A combination of identifiers + coordinates is allowed, but only one unique value per source is needed.Note: script resolves names with Simbad by default. Future change should be to allow the user to utilise the "name_resolver" field in the Perl script, assuming this is possible in the Java one?

Coordinates can be given in the following forms:
- DEGREES: e.g., "105.0 54.7"
- SEXAGESIMAL: e.g., "12 29 06.70,02 03 08.6"
Make sure the coordinates are one per line in the form "RA DEC".

- Method:
This script parallelizes the search list if a long list is given. The source list is split into chunks with a maximum number of 10 sources per chunk. The script then parallelizes queries over each chunk and stitches the results back together at the end.

- User-defined tables:
Any tables available on HEASARC can be queried. The user just needs to add  these to the "INSTRUMENTS.csv" file present in the folder, along with the  instrument, instrument offset (in arcmin) and name of the exposure column (can be found by clicking on the corresponding instrument and table here:  https://heasarc.gsfc.nasa.gov/cgi-bin/W3Browse/w3browse.pl). Note: the constraints column can be left blank, or the user can add their  own (at their own risk!).

- Summary file:
The script finishes by creating a summary file -- "SUMMARY.csv" that gives the details of all X-ray data found for every source, in the original order provided by the user.


- to make this script work on other machines:
Change setupDir to the directory containing users.jar and INSTRUMENTS.csv Then change to the directory containing your source list and use the command: >>> xq -t my_sources.txt -w y (-w y makes sure wget commands are found) where my_sources.txt is the file containing the source list

MAKE SURE coordinate RA and dec are separated by a space!

- other links:
Full user guide for the Java script: https://heasarc.gsfc.nasa.gov/xamin/doc/CLIUsersGuide.html
