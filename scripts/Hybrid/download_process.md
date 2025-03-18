# download.sh documentation
## Process to download scRNAseq
<ol>
    <li>Change the execute rights of the bash script </br> 
    <code>chmod +x ./download.sh </code></li>
    <li>Modify the csv file with the urls, hash types and checksums necessary for your download</li> 
    <li>Run script </br> 
    <code>./download.sh </code> </li>
    <li>Find your files in data directory</li> 
    
</ol>


## download.sh pseudocode
**INPUT:** urls.csv </br>
**OUTPUT:** downloaded files </br> 
The script <i>download.sh</i> does the following steps: </br>
<ol>
    <li>Obtain url, hash type and checksum from urls.csv </li>
    <li>Run in the background a function named download_and_extract that does the following</li> 
    <ol>
        <li>Download a file using url</li>
        <li>Check the checksum depending on the hash type</li>
        <li>Extract the files depending on the type (tar, tar.gz, tar.bz2, tar.xz, zip)</li>
        <li>Remove the tar/zip files</li>
    </ol>   
<ol>