##### Making a database - mcrA edition #####
# Author: Stephania L. Tsola
# with help from Dr. Kate Randall and Dr. Yizhu Zhu
# Date: Feb 2022



 ## Download reference database from fungenes
# Find mcrA gene -> "Begin Analysis" -> nucleotide download -> download labels with accession number -> fasta format -> unaligned seqs -> Doqnload
##Keep in mind that you can only download 10000 sequences from fungenes at a time so you migh have to download them in chunks. To combine all files into one put them in a seperate folder and:
cat * > mynewfile.fasta  #the * means "take all files in the folder"

## Remove header from fungene fasta file and put in other file
grep -e ">" fungene_all.fasta > fungene_all_Headers.txt

## e.g >CP008746  location=complement(3239792..3241504),organism=Methanosarcina barkeri CM1,definition=methyl-coenzyme M reductase alpha subunit McrA

#I want to remove everything besides the first code "CP008746"
python3
with open("fungene_all_Headers.txt","r+") as data_file:
    list_of_elements=[] # creating an empty list to store all elements
    for line in data_file:
        data=line.split()
        del data[1:]
        print(str(data))
        list_of_elements.append(str(data)) # appending the list with the elements
#above basically converts all lines into lists (1st elements of each line one list, 2nd element other list etc) then I tell it to delete all lists besides the first (0 in python) which contain what I want
with open('fungene_AllAccNum.txt', 'w') as f: # writing to a file
    for line in list_of_elements:
        f.write(line)
        f.write('\n')
# Because I made lists when everything gets put into a new file you end up with "['>FR871724']" but with find and replace in the new text file you can get rid of the extra spaces, [], > and ' ---> Fix code to do this automatically

# Splitting One Big text File Into Multiple Files - we need the accesion numbers to query ncbi but ncbi is very strict with number of queries so need to split the file for it to work
split -l 200 mybigfile.txt

#Add a header to each file "#SampleID" otherwise next steps won't work (qiime won't recognise them)

# Loop to parse through all the individual files and query ncbi for the accession numbers - Use RESCRIpt plugin for qiime
  #Put all individual files in their own folder and run script within the folder
  #RESCRIPt: https://library.qiime2.org/plugins/rescript/27/

conda activate qiime2-2021.11

for file in *; do qiime rescript get-ncbi-data --m-accession-ids-file ${file} --o-sequences seqs/${file}-seqs.qza --o-taxonomy tax/${file}-tax.qza; done

#Put all seqs files and tax files in their own folders and once within run code below to merge all the individual tax and seqs files into one big seqs and tax file
qiime feature-table merge-taxa --i-data *-tax.qza --o-merged-data all-tax.qza

qiime feature-table merge-seqs --i-data *-seqs.qza --o-merged-data All-seqs.qza

# Unzip the big seqs and tax files to check everything is ok with the data
unzip all-tax.qza
unzip All-seqs.qza

### Filtering of tax and seqs files. After making a note of what you want to filter out you can use qiime filter-seqs to filter out stuff from the seqs files. The tax file you need to do "manually" from the text file. At the end make sure the two files match!! I wonder what happens if they don't match - I guess the classifier training is just matching the two files together. Do the extra tax get discarded or does qiime have an error?

###Removing from seqs file to match the tax file ### Can I put all these together by running the p-exclude flag multiple times?
qiime taxa filter-seqs \
--i-sequences All-seqs.qza \
--i-taxonomy all-tax.qza \
--p-exclude k__Bacteria \
--o-filtered-sequences AllnoBac.qza

qiime taxa filter-seqs \
--i-sequences AllnoBac.qza \
--i-taxonomy all-tax.qza \
--p-exclude s__archaeon \
--o-filtered-sequences NoBacUncArc.qza

qiime taxa filter-seqs \
--i-sequences NoBacUncArc.qza \
--i-taxonomy all-tax.qza \
--p-exclude f__Archaea \
--o-filtered-sequences NoBacUncArc2.qza

qiime taxa filter-seqs \
--i-sequences NoBacUncArc2.qza \
--i-taxonomy all-tax.qza \
--p-exclude f__Euryarchaeota \
--o-filtered-sequences NoBacUncArc3.qza

qiime taxa filter-seqs \
--i-sequences NoBacUncArc3.qza \
--i-taxonomy all-tax.qza \
--p-exclude k__Environmental \
--o-filtered-sequences NoBacUncArcEnv.qza

qiime taxa filter-seqs \
--i-sequences NoBacUncArcEnv.qza \
--i-taxonomy all-tax.qza \
--p-exclude g__archaeon \
--o-filtered-sequences NoBacUncArc4Env.qza

qiime taxa filter-seqs \
--i-sequences NoBacUncArc4Env.qza \
--i-taxonomy all-tax.qza \
--p-exclude g__methanogenic \
--o-filtered-sequences NoBacUncArc5Env.qza

###Removing from taxonomy file! ### after unzipping qza file and converting file back to txt from tsv (renaming works)
# Remove bacteria
sed '/k__Bacteria/d' All_NCBI_tax.txt > NoBac.txt
# Remove uncultured archaea1
sed '/s__archaeon/d' NoBac.txt > NoBacUncArc.txt
# Remove uncultured archaea2
sed '/f__Archaea/d' NoBacUncArc.txt > NoBacUncArc2.txt
# Remove uncultured archaea3
sed '/f__Euryarchaeota/d' NoBacUncArc2.txt > NoBacUncArc3.txt
# Remove Environmental samples
sed '/k__Environmental/d' NoBacUncArc3.txt > NoBacUncArcEnv.txt
# Remove uncultured archaea4
sed '/g__archaeon/d' NoBacUncArcEnv.txt > NoBacUncArc4Env.txt
# Remove uncultured archaea5
sed '/g__methanogenic/d' NoBacUncArc4Env.txt > NoBacUncArc5Env.txt

#Remove header from final txt file otherwise next step won't work (there must be code that goes around this but couldn't find it plus it was easy to just remove the header by deleting manually)

#Reimport tax file
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path NoBacUncArc5Env-tax.txt \
--output-path NoBacUncArc5Env-tax.qza

# “Culling” low-quality sequences with cull-seqs with rescript
#Removes sequences that contain 5 or more ambiguous bases (IUPAC compliant ambiguity bases) and any homopolymers that are 8 or more bases in length
qiime rescript cull-seqs \
    --i-sequences NoBacUncArc5Env-seqs.qza \
    --o-clean-sequences NoBacUncArc5Env-seqs-cleaned.qza
    
#Trains classifier
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  NoBacUncArc5Env-seqs-cleaned.qza \
  --i-reference-taxonomy NoBacUncArc5Env-tax.qza \
  --o-classifier NoBacUncArc5Env-tax-classifier.qza

#Everything is now ready to use in qiime2

