#Code for download data
dir.create("1_Get_Data", showWarnings = FALSE)
dir.create("post", showWarnings = FALSE)
dir.create("pre", showWarnings = FALSE)

# Get pre Data
download.file("https://cf.10xgenomics.com/samples/cell-exp/1.1.0/aml035_pre_transplant/aml035_pre_transplant_fastqs.tar",
              destfile = "pre/aml035_pre_transplant_fastqs.tar", method = "wget")
# Untar
untar("pre/aml035_pre_transplant_fastqs.tar", exdir = "pre")

# Get post data
download.file("https://cf.10xgenomics.com/samples/cell-exp/1.1.0/aml035_post_transplant/aml035_post_transplant_fastqs.tar",
              destfile = "post/aml035_post_transplant_fastqs.tar", method = "wget")

# untar
untar("post/aml035_post_transplant_fastqs.tar", exdir = "post")
