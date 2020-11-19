"""
Author: Peter Boorman, Czech Academy of Sciences
This script searches the High Energy Astrophysics Science Archive Research Center 
(HEASARC) archive for any archival X-ray data available for a user-provided list 
of targets.

- Input:
Targets can be provided in a text file, with one source identifier (or set of 
coordinates) per line. A combination of identifiers + coordinates is allowed, 
but only one unique value per source is needed.
Note: script resolves names with Simbad by default. Future change should be to
allow the user to utilise the "name_resolver" field in the Perl script, 
assuming this is possible in the Java one?

Coordinates can be given in the following forms:
- DEGREES: e.g., "105.0 54.7"
- SEXAGESIMAL: e.g., "12 29 06.70,02 03 08.6"
Make sure the coordinates are one per line in the form "RA DEC".

- Method:
This script parallelizes the search to avoid overwhelming the servers used 
to search. Thus currently, the source list is split into chunks with a 
maximum number of 10 sources per chunk. The script then parallelizes queries 
over each chunk and stitches the results back together at the end.

- User-defined tables:
Any tables available on HEASARC can be queried. The user just needs to add 
these to the "INSTRUMENTS.csv" file present in the folder, along with the 
instrument, instrument offset (in arcmin) and name of the exposure column
(can be found by clicking on the corresponding instrument and table here: 
https://heasarc.gsfc.nasa.gov/cgi-bin/W3Browse/w3browse.pl).
Note: the constraints column can be left blank, or the user can add their 
own (at their own risk!).

- Summary file:
The script finishes by creating a summary file -- "SUMMARY.csv" that 
gives the details of all X-ray data found for every source, in the 
original order provided by the user.


- to make this script work on other machines:
Change setupDir to the directory containing users.jar and INSTRUMENTS.csv
Then change to the directory containing your source list and use the command:
>>> xq -t my_sources.txt -w y
(-w y makes sure wget commands are found)
where my_sources.txt is the file containing the source list

MAKE SURE coordinate RA and dec are separated by a space!

- other links:
Full user guide for the Java script: https://heasarc.gsfc.nasa.gov/xamin/doc/CLIUsersGuide.html
"""

import os,sys,argparse,datetime,glob,multiprocessing,uuid,re,requests
import pandas as pd
import numpy as np
from tqdm import *
from astropy.time import Time

setupDir = "/Users/pboorman/Dropbox/data/0_analysis/xrayQuery/"
workingDir = os.getcwd()

def createSummary(outputDir,srcList,TABLES):
    """
    Create summary file listing the data found for all sources in the original order
    they were given in. Sources with no data found are included in the summary for
    clarity
    """
    print("Creating summary file for all sources in original order given...")
    fullSrcList=pd.read_csv("../"+srcList,names=["0INPUT"])
    fullSrcList.loc[:,"1MATCHED_TARGET"]=np.repeat("NOT_FOUND",len(fullSrcList))
    ins_abbreviations=["NUS","CHA","SUZ","XRT","XMM"]
    totObs=0.
    totTime=0.
    with tqdm(total=len(TABLES)*len(fullSrcList)) as pbar:
        for i,table in enumerate(TABLES):
            INSTR=ins_abbreviations[i]
            nObs="N_%(INSTR)s" %locals()
            medOffset="medSep_%(INSTR)s" %locals()
            medExp="medks_%(INSTR)s" %locals()
            totExp="totks_%(INSTR)s" %locals()
            fullSrcList.loc[:,nObs]=int(0)
            fullSrcList.loc[:,medOffset]=0.
            fullSrcList.loc[:,medExp]=0.
            fullSrcList.loc[:,totExp]=0.
            for s,src in fullSrcList.iterrows():
                tempInst=pd.read_csv("%(table)s.csv" %locals(), dtype = "str")
                specificObs=tempInst.loc[tempInst["b_target"]==src["0INPUT"]]
                if len(specificObs)>0:
                    fullSrcList.loc[s,"1MATCHED_TARGET"]=specificObs["a_name"].values[0]
                fullSrcList.loc[s,nObs]=str(len(specificObs))
                fullSrcList.loc[s,medOffset]=np.round(specificObs["OFFSET_ARCMIN"].astype(float).median(),2)
                fullSrcList.loc[s,medExp]=np.round(specificObs["EXPOSURE"].astype(float).median()/10**3,2)
                fullSrcList.loc[s,totExp]=np.round(specificObs["EXPOSURE"].astype(float).sum()/10**3,2)
                pbar.update()
            totObs+=fullSrcList[nObs].astype(int).sum()
            totTime+=fullSrcList[totExp].sum()
            
    avObs=np.round(totObs/len(fullSrcList),1)    
    fullSrcList=fullSrcList.reindex(sorted(fullSrcList.columns), axis=1)
    fullSrcList=fullSrcList.rename(columns={
        "0INPUT":'INPUT',
        "1MATCHED_TARGET":'MATCHED_TARGET'
        })
    fullSrcList.to_csv("SUMMARY.csv" %locals(),na_rep="NaN",index=False)
    print("Done.")
    return int(totObs),avObs,np.round(totTime,2)

def getInstrument_df(chunk_df):
    """
    This generates a dataframe from the tables specified in the external
    "INSTRUMENTS.csv" file.
    Users can add additional tables here to search over.
    Additional values are appended to the dataframe to use later.
    """
    ID=str(uuid.uuid4().hex)
    info_df=pd.read_csv(setupDir + "/INSTRUMENTS.csv").fillna("")
    info_df.loc[:,"constraint"]=""
    info_df.loc[info_df.instrument=="nustar","constraint"]="constraint='exposure_a > 1000'"
    info_df=info_df.loc[info_df.table==chunk_df.TABLE.values[0]]
    javaPath = setupDir + "/users.jar"
    javaCmd="java -jar %(javaPath)s" %locals()
    cwd=os.getcwd()
    info_df.loc[:,"javaUpload"]=cwd+"/"+ID+"_sources.txt"
    info_df.loc[:,"javaOutput"]=cwd+"/"+ID+"_"+info_df["table"]+".xls"
    
    info_df.loc[:,"JAVA_CMD"]=javaCmd + \
    " table="+info_df["table"].astype(str) + \
    " upload="+info_df["javaUpload"].values[0] + \
    " offset=a:b:"+info_df["offset_arcmin"].astype(str) + \
    " showoffsets " + \
    info_df["constraint"].astype(str) + \
    " fields=" + columns + " " + \
    "output="+info_df["javaOutput"].astype(str) + \
    " format=excel"
    chunk_df["0INPUT"].to_csv(info_df["javaUpload"].values[0],header=False,index=False)
    ## alternatively, fields=standard will return other columns
    return info_df

def chunkSearch(chunk_df):
    """
    This searches for data on a specific chunk of the original source list.
    First save source list chunk to unique name.
    ... Then create instrument df with java command featuring the specific unique name.
    """

    instrument_df=getInstrument_df(chunk_df)
    runJavaCmd=instrument_df.JAVA_CMD.values[0]
    javaUpload=instrument_df.javaUpload.values[0]
    javaOutput=instrument_df.javaOutput.values[0]
    instrument=instrument_df.instrument.values[0]
    os.system(runJavaCmd)
    xlsResFile=instrument_df.javaOutput.values[0]
    if os.stat(xlsResFile).st_size != 0:
        temp=pd.read_excel(xlsResFile,skiprows=1, dtype = "str").dropna(subset=["b_target"])
        temp=temp.rename(columns={
            instrument_df.exposure_column.values[0]:'EXPOSURE',
            "a_obsid":'OBSID',
            "_delta_ab":'OFFSET_ARCMIN',
            "_delta_ba":'OFFSET_ARCMIN'})
        csvResFile=xlsResFile.replace(".xls",".csv")
        temp.to_csv(csvResFile,index=False)
        os.remove(xlsResFile)
        


def parallelizeQueries(searchNames_df,FUNCTION=chunkSearch,NCORES=4):
    """
    Split entire source list into chunks of maximum length maxChunkLen.
    The function then parallelizes the queries of multiple chunks
    depending on how many cores are available.
    """
    maxChunkLen=10 ## maximum number of rows in a dataframe chunk
    nChunk=int(np.ceil(len(searchNames_df)/maxChunkLen))
    tab=searchNames_df.TABLE.values[0]
    print("Searching %(tab)s..." %locals())
    DATA_SPLIT=np.array_split(searchNames_df,nChunk)
    with multiprocessing.Pool(processes=NCORES) as p:
        with tqdm(total=nChunk) as pbar:
            for i,_ in enumerate(p.imap_unordered(FUNCTION,DATA_SPLIT)):
                pbar.update()
    for file in glob.glob("*_sources.txt"):
        os.remove(file)    

def create_obstime_column(tab, df_in):
    """
    Used to convert the MJD date of observation to human-readable form
    """
    if tab == "swiftmastr":
        time_column = "a_start_time"
    else:
        time_column = "a_time"

    df = df_in.copy()
    df.loc[:, "obs_year"] = df.dropna(subset = [time_column])[time_column].apply(lambda x: Time(float(x), format = "mjd").to_datetime().year)
    df.loc[:, "obs_month"] = df.dropna(subset = [time_column])[time_column].apply(lambda x: Time(float(x), format = "mjd").to_datetime().month)
    df.loc[:, "obs_day"] = df.dropna(subset = [time_column])[time_column].apply(lambda x: Time(float(x), format = "mjd").to_datetime().day)

    df.loc[:, "OBSTIME"] = df.dropna(subset = [time_column])[time_column].apply(lambda x: Time(float(x), format = "mjd").to_datetime().strftime("%Y%b%d"))
    return df

def startQueries(searchFile, wgetLocs = "n", filterData = "y"):
    """
    Main function for running the query.
    Reads the input target list, and starts the parallelization.
    Finishes by creating summary file for every target searched,
    irrespective of whether any data was found or not.
    """
    start=datetime.datetime.now()    
    instruments=pd.read_csv(setupDir + "/INSTRUMENTS.csv")
    ## uncomment to only find data for a select instrument
    # instruments = instruments.loc[instruments["instrument"] == "xmm"]
    tables=instruments.table

    for i,tab in enumerate(tables):
        searchNames_df=pd.read_csv("../" + searchFile,names=["0INPUT"])
        searchNames_df.loc[:,"TABLE"]=tab
        
        nCores=multiprocessing.cpu_count()
        if nCores>4:nCores-=4
        parallelizeQueries(searchNames_df,chunkSearch,nCores)
    
        df=pd.DataFrame()
        for resFile in glob.glob("*%(tab)s*" %locals()):
            if os.path.getsize(resFile):
                temp=pd.read_csv(resFile, dtype = "str")
                df=pd.concat([df,temp],ignore_index=True,sort=True)
            os.remove(resFile)
        df = create_obstime_column(tab, df)
        df.to_csv("%(tab)s.csv" %locals(),index=False)
        if filterData == "y":
            ## note Suzaku and Swift/XRT do not have the a_status column
            ## (assuming this is because all data is available).
            if "a_status" in df.columns:
                df = df.loc[df["a_status"].astype(str) == "archived"]
                df.to_csv("%(tab)s.csv" %locals(),index=False)
        if wgetLocs=="y":getFileLocFTP(tab)
        print("Done.")

    Nobs,avObs,Ntime=createSummary(outputDir,searchFile,tables)
    stop=datetime.datetime.now()
    timeTaken=np.round((stop-start).total_seconds()/60.,1)
    Nsrcs=len(searchNames_df)
    print("\nA total of %(Nobs)s observations (average %(avObs)s per source) \
covering %(Ntime)s ks were found in %(timeTaken)s mins for %(Nsrcs)s sources." %locals())

    
def exists(path):
    """
    Checks if data is available at FTP location
    """
    r=requests.head(path)
    return r.status_code==requests.codes.ok


def getFileLocFTP(table):
    """
    Very useful links:
    https://heasarc.gsfc.nasa.gov/docs/FTPWarning.html

    """
    dataPaths={
    "numaster":"nustar/data/obs/OBSID[1:3]/OBSID[0]/OBSID/",
    "chanmaster":"chandra/data/byobsid/OBSID[-1]/OBSID/",
    "suzamaster":"suzaku/data/obs/OBSID[0]/OBSID/",
    "swiftmastr":"/swift/data/obs/YEAR_MONTH.zfill(2)/OBSID/xrt/",
    "xmmmaster":"xmm/data/rev0/OBSID/ODF/"}

    ## number of dirs to put in wget command
    ## these ensure that the data will be downloaded with an
    ## encompassing obsid folder
    ndirs = {
    "numaster": "6",
    "chanmaster": "5",
    "suzamaster":"5",
    "swiftmastr":"5",
    "xmmmaster":"4"}
    
    info_df=pd.read_csv(setupDir + "/INSTRUMENTS.csv").fillna("")
    inst_df=info_df.loc[info_df.table==table]
    data=pd.read_csv("%(table)s.csv" %locals(), dtype = "str")
    
    if len(data)==0: return
    # if table=="swiftmastr":return

    data.loc[:,"dataLoc"]=""
    data.loc[:,"wget_cmd"]=""
    wget_start = "wget\
 -nH\
 --no-check-certificate\
 --cut-dirs=X\
 -r\
 -l0\
 -c\
 -N\
 -np\
 -R\
 'index*'\
 -erobots=off\
 --retr-symlinks "
    print("Finding FTP locations for %(table)s.csv" %locals())
    with tqdm(total=len(data)) as pbar:
        for i,obs in data.iterrows():
            try:
                dataPath="https://heasarc.gsfc.nasa.gov/FTP/"+dataPaths[table]
            except:
                print("Unknown file locations for %(table)s. Skipping..." %locals())
                return
            OBSID=str(obs.OBSID)
            tryPath=dataPath.replace("OBSID[0]",OBSID[0])
            tryPath=tryPath.replace("OBSID[1:3]",OBSID[1:3])
            tryPath=tryPath.replace("OBSID[-1]",OBSID[-1])
            # tryPath=tryPath.replace("OBSID.zfill(10)",OBSID.zfill(10)) ## don't need to zfill now that obsids are read in as strings
            tryPath=tryPath.replace("YEAR_MONTH.zfill(2)", str(obs["obs_year"]) + "_" + str(obs["obs_month"]).zfill(2))
            tryPath=tryPath.replace("OBSID",OBSID)
            tryPath=re.sub('\/+', '/',tryPath).replace(":/","://")
            if exists(tryPath):
                fullDataLoc=tryPath
                full_wget_cmd = wget_start.replace("cut-dirs=X", "cut-dirs=" + ndirs[table])
            else:
                fullDataLoc="UNKNOWN"
                full_wget_cmd = "UNKNOWN"

            data.loc[i,"dataLoc"]=fullDataLoc
            data.loc[i,"wget_cmd"]=full_wget_cmd + " " + fullDataLoc
            pbar.update()  
    # pd.options.display.max_colwidth=100
    data.to_csv("%(table)s.csv" %locals(),index=False)



parser=argparse.ArgumentParser(
    description="Find archival X-ray observations from a list of source \
names and/or coordinates.")

parser.add_argument('--textfile',
    '-t',
    nargs='?',
    const="mysources.txt",
    default="mysources.txt",
    type=str,
    help="Name of textfile containing names and/or coordinates, one \
per line. Default textfile name == 'mysources.txt'")

parser.add_argument('--wgetlocs',
    '-w',
    nargs='?',
    const="n",
    default="n",
    type=str,
    help="Choice to find the wget locations of all data found. Warning: \
takes a long time. Can be y|Y|n|N, default == n.")

parser.add_argument('--filterdata',
    '-f',
    nargs='?',
    const="n",
    default="n",
    type=str,
    help="Choice to filter the observations found for data that is archived, \
i.e. available. Can be y|Y|n|N, default == n.")

parser.add_argument('--columns',
    '-c',
    nargs='?',
    const="all",
    default="all",
    type=str,
    help="Choice to include 'standard' columns or 'all' columns in returned tables.")


args=parser.parse_args()
searchFile=args.textfile
wgetLocs=args.wgetlocs
filterData = args.filterdata
columns = args.columns

if wgetLocs.lower()!="y":wgetLocs="n"
if filterData.lower()!="n":wgetLocs="y"
if (columns.lower() != "all") and (columns.lower() != "standard"): columns = "all"
if not os.path.isfile(searchFile):
    searchFile = workingDir + "/" + searchFile

if not os.path.isfile(searchFile):
    print("File %(searchFile)s does not exist. Exiting..." %locals());sys.exit()

outputDir=datetime.datetime.now().strftime("xq_%Y%b%d_%H%M%S")
if not os.path.isdir(outputDir):
    os.mkdir(outputDir)
else:
    print("%(outputDir)s already exists. Exiting..." %locals());sys.exit()
os.chdir(outputDir)
startQueries(searchFile,wgetLocs, filterData)
os.chdir("../")

