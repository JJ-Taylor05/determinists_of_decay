# Script to test that CMfinder is functioning 
# Run this in the same directory as the xrrna_flavivirus_3utr.fa file (available in my Github repo "JJ-Taylor05/determinists_of_decay/test_files/"

mkdir cmfinder_test

cmfinder.pl xrrna_flavivirus_3utr.fa

mv *.motif.* *.cm.* cmfinder_test/
