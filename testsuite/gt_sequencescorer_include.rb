Name "gt_sequencescorer_check_gramdistance"
Keywords "gt_sequencescorer qgram distance"
Test do
  run "#{$bin}gt scorer -qgram -s #{$testdata}EncodeTransProt7_One.txt #{$testdata}EncodeTransProt7_Two.txt -q 4"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt7_q4.out"
  run "#{$bin}gt scorer -qgram -s #{$testdata}EncodeTransProt11_One.txt #{$testdata}EncodeTransProt11_Two.txt -q 5"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt11_q5.out"
end 

Name "gt_sequencescorer_check_fscore"
Keywords "gt_sequencescorer fscore distance"
Test do
  run "#{$bin}gt scorer -fscore -s #{$testdata}EncodeTransProt7_One.txt #{$testdata}EncodeTransProt7_Two.txt -k 4"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt7_k4.out"
  run "#{$bin}gt scorer -fscore -s #{$testdata}EncodeTransProt11_One.txt #{$testdata}EncodeTransProt11_Two.txt -k 5"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt11_k5.out"
end 
