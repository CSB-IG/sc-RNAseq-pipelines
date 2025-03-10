# download_script.py documentation

## Process to download scRNAseq
<ol>
    <li>Import conda environment</br> 
    <code>conda env create -f environment.yml</code></li>
    <li>Modify the csv file with the urls, hash types and checksums necessary for your download</br> 
    <code>conda activate hackaton-sc</code></li>
    <li>Run script </br> 
    <code>python download_scripts.py </code> </li>
    <li>Find your files in data directory</li> 
    
</ol>


## download_script.py pseudocode
The pseudocode is similar as with Hybrid script. 
**INPUT:** urls.csv </br>
**OUTPUT:** downloaded files </br> 
The script <i>download_script.py</i> does the following steps: </br>
<ol>
    <li>Obtain url, hash type and checksum from urls.csv </li>
    <li>Run a function named download_and_verify that does the following</li> 
    <ol>
        <li>Download a file using url (raise exceptions with http responses)</li>
        <li>Check the checksum depending on the hash type with function verify_checksum</li>
        <li>Extract the files depending on the type (tar, tar.gz, tar.bz2, tar.xz, zip)</li>
        <li>Remove the tar/zip files</li>
    </ol>   
<ol>

## Functions
### calculate_checksum()
**DESCRIPTION** Function to calculate checksum of a file with 3 hash types [SHA256, SHA1, MD5], to later check if file's integrity is not compromised.
**INPUT:** </br>
<ul>
    <li><i>file_path (string)</i>: Path of the directory in which the downloaded file is located</li>
    <li><i>hash_type (string [SHA256, SHA1, MD5])</i>: Type of hashing needed to check file's integrity</li>
</ul> 
**OUTPUT** </br>
<ul>
    <li><i>hash_func().hexdigest (string)</i>: String of checksum to compare with website's checksum and verify if file is intact.</li>
    <li><i>None</i>: It returns None when the hash_type is not supported, you can modify the function to include the desired hash_type</li>
</ul> 

### verify_checksum()
**DESCRIPTION** Function that uses the function calculate_checksum to verify if file's integrity is intact. In case that the file is corrupted, the file will be deleted.
**INPUT:** </br>
<ul>
    <li><i>file_path (string)</i>: Path of the directory in which the downloaded file is located</li>
    <li><i>hash_type (string [SHA256, SHA1, MD5])</i>: Type of hashing needed to check file's integrity</li>
    <li><i>expected_checksum (string)</i>: Set of characters obtained from the website to check for file's integrity (need to input checksum and hash in urls.csv)</li>

</ul> 
**OUTPUT** </br>
<ul>
    <li><i>bool</i>: True if checksum match, false if it doesn't match.</li>
</ul> 

### extract_file()
**DESCRIPTION** Function that extracts files deepending of its type (tar, tar.gz, tar.bz2, tar.xz, zip)
**INPUT:** </br>
<ul>
    <li><i>file_path (string)</i>: Path of the directory in which the downloaded file is located</li>
    <li><i>dir_path (string)</i>: Directory in which all the extracted files will be saved</li>

</ul> 
**OUTPUT** </br>
<ul>
    <li><i>None</i></li>
</ul> 


### download_and_verify()
**DESCRIPTION** Function that extracts files deepending of its type (tar, tar.gz, tar.bz2, tar.xz, zip)
**INPUT:** </br>
<ul>
    <li><i>file_path (string)</i>: Path of the directory in which the downloaded file is located</li>
    <li><i>dir_path (string)</i>: Directory in which all the extracted files will be saved</li>

</ul> 
**OUTPUT** </br>
<ul>
    <li><i>None</i></li>
</ul> 
