source("tcr_plots.r")

# In this example we have three files: 0XX1A_P, 0XX1A_B30, 0XX1A_B123.
# These files represent different time-points of a specific patient.
# These time-points have been encoded as "sampletype".
# We generate a list of files, "lfiles"


patient_id="0XX1"
chain="A"
title=paste0("Patient: ",patient_id," chain: ",chain)

# the graphs will appear with the same order of the sampletype
sampletype=c("P", "B30", "B123")

# the list of files
lfiles=paste(paste0(patient_id, chain), sampletype, sep="_")

basicnf="_sample_0XX1_"

## databind contains the informations of all the patients
## its columns are:  trxvseqtrxj, frequency, patient
## these columns represent, repsectively, the HTRXV_CDR3_TRXJ sequence,
## the frequency of such sequence, the patient from which that sequence comes from
databind=create_df(lfiles, sampletype)


## we need a list of dataframe, one for each patient, although this is a redundant with databind
ldataf=vector("list", length(lfiles))
for(i in 1:length(sampletype)){
  ldataf[[i]]= databind[databind$patient==sampletype[i], ]
}


# one or more sequence that we want to follows in the various timepoint
cseq=c("hTRAV27_CAGARGGADGLTFG_hTRAJ45", "hTRAV17_CATDQAGTALIFG_hTRAJ15")


## trace_seq will generate a dataframe for each patient/timepoint for a specific sequence.
## If a sequence is found, then it associates its frequency at that given timepoint
## otherwise the frequency is zero
tddf=trace_seq(ldataf, cseq, sampletype)

# create the violin plot for our dataframe
bv=base_violin(databind, title)

# we trace the sequences we wanted (remember, for scale color manual we need 
# as many colors as sequences we want to trace)
p=bv+geom_line(data=tddf, aes(x=patient, y=frequency, group=trxvseqtrxj, colour=trxvseqtrxj))+
  scale_color_manual(values=c("green", "blue"), name="TCR sequence")

ggsave(plot=p, filename=paste0("graph_", basicnf,".jpg"), device="jpg", width=14)

## we want to trace the first 10 (and 20) sequences of sampletype[1]( that corresponds to ldataf[[1]])
## that are also present in the other samples.
## if a sequence has frequency ZERO we will ignore it and check if the next sequence is present in all the samples
nseq=c(10, 20)
for(n in nseq){
    lseq=find_n_common_seq(ldataf, n, ldataf[[1]])

    tddf2=trace_seq(ldataf, lseq, sampletype)
    
    r=bv+geom_line(data=tddf2, aes(x=patient, y=frequency, group=trxvseqtrxj), color="black",
                   size=0.25, alpha=0.5)
    
    targetnf=paste0("first_",n,"_of_", sampletype[1],"_")
    
    nf=paste0(targetnf, basicnf, ".jpg")
    ggsave(plot=r, filename=nf, device="jpg", width=14)

}
