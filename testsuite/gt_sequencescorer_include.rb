Name "gt_sequencescorer_check_gramdistance_TransProt7"
Keywords "gt_sequencescorer scorer qgram TransProt7"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerOne.txt"
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerTwo.txt"
  run "#{$bin}gt scorer -qgram -s SequencescorerOne.txt SequencescorerTwo.txt -q 4"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt7_q4.out"
end

Name "gt_sequencescorer_check_gramdistance_TransProt11"
Keywords "gt_sequencescorer scorer qgram TransProt11"
Test do
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}SequencescorerOne.txt"
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}SequencescorerTwo.txt"
  run "#{$bin}gt scorer -qgram -s SequencescorerOne.txt SequencescorerTwo.txt -q 5"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt11_q5.out"
end

Name "gt_sequencescorer_check_fscore_TransProt7"
Keywords "gt_sequencescorer scorer fscore TransProt7"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerOne.txt"
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerTwo.txt"
  run "#{$bin}gt scorer -fscore -s SequencescorerOne.txt SequencescorerTwo.txt -k 4"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt7_k4.out"
  end
  
Name "gt_sequencescorer_check_fscore_TransProt11"
Keywords "gt_sequencescorer scorer fscore TransProt11"
Test do
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}SequencescorerOne.txt"
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}SequencescorerTwo.txt"  
  run "#{$bin}gt scorer -fscore -s SequencescorerOne.txt SequencescorerTwo.txt -k 5"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt11_k5.out"
end 


Name "gt_sequencescorer_check_gramdistance_lowerbound"
Keywords "gt_sequencescorer scorer qgram lowerbound"
Test do
  outfile = File.new("output.txt", "w")
  infileone = File.new("#{$testdata}SequencescorerOne.txt", "r")
  infiletwo = File.new("#{$testdata}SequencescorerTwo.txt", "r")
  onearray = infileone.readlines
  twoarray = infiletwo.readlines
  infileone.close
  infiletwo.close
  for i in 0...onearray.length
    if(onearray[i][">"])
      next
    end
    for j in 0...twoarray.length
      if(twoarray[j][">"])
        next
      end
      run "#{$bin}gt dev paircmp -e -ss #{onearray[i].chomp} #{twoarray[j].chomp}"
      outfile << last_stdout << "\n"
    end
  end
  outfile.close
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerOne.txt"
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerTwo.txt"
  run "#{$bin}gt scorer -qgram -s SequencescorerOne.txt SequencescorerTwo.txt -q 4"
  outfile = File.new("#{last_stdout}", "r")
  qgramarray = outfile.readlines
  outfile.close
  for i in 0...qgramarray.length
    if(qgramarray[i][">"])
        next
    end
    tmp = (qgramarray[i].split)
    qgramarray.delete(qgramarray[i])
    qgramarray.insert(i,tmp[tmp.length-1])
  end
  outfile = File.new("output.txt", "r")
  edistarray = outfile.readlines
  outfile.close
  for i in 0...edistarray.length
    tmp = File.new("#{edistarray[i].chomp}", "r")
    tmparray = tmp.read.split()
    if("#{tmparray[0]}".to_i < ("#{qgramarray[i]}".to_i)/8)
      exit 1
    end
  end
end 
