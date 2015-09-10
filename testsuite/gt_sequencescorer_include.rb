Name "gt sequencescorer Alphabet"
Keywords "gt_sequencescorer scorer alphabet"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerOne.txt"
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}SequencescorerTwo.txt"
  run "#{$bin}gt scorer -qgram -ii SequencescorerOne.txt "\
      " SequencescorerTwo.txt -q 4", :retval => 1
end

Name "gt sequencescorer Fileformat"
Keywords "gt_sequencescorer scorer Fileformat"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerOne.txt"
  run "#{$bin}gt scorer -qgram -ii SequencescorerOne.txt "\
      "SequencescorerTwo.txt -q 4", :retval => 1
  run "#{$bin}gt scorer -qgram -ii -q 4", :retval => 1
  run "#{$bin}gt scorer -qgram -ii SequencescorerOne.txt "\
      "SequencescorerOne.txt SequencescorerOne.txt -q 4", :retval => 1
end

Name "gt sequencescorer Qgramdistance"
Keywords "gt_sequencescorer scorer qgram"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerOne.txt"
  run "#{$bin}gt scorer -qgram -ii SequencescorerOne.txt -q 4", :retval => 0
end

Name "gt sequencescorer Qgramdistance TransProt7"
Keywords "gt_sequencescorer scorer qgram TransProt7"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerOne.txt"
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerTwo.txt"
  run "#{$bin}gt scorer -qgram -ii SequencescorerOne.txt "\
      "SequencescorerTwo.txt -q 4"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt7_q4.out"
end

Name "gt sequencescorer Qgramdistance TransProt11"
Keywords "gt_sequencescorer scorer qgram TransProt11"
Test do
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}SequencescorerOne.txt"
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}SequencescorerTwo.txt"
  run "#{$bin}gt scorer -qgram -ii SequencescorerOne.txt "\
      "SequencescorerTwo.txt -q 5"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt11_q5.out"
end

Name "gt sequencescorer Fscore"
Keywords "gt_sequencescorer scorer fscore"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerOne.txt"
  run_test "#{$bin}gt scorer -fscore -ii SequencescorerOne.txt -k 6", :retval => 0
end

Name "gt sequencescorer Fscore TransProt7"
Keywords "gt_sequencescorer scorer fscore TransProt7"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerOne.txt"
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerTwo.txt"
  run "#{$bin}gt scorer -fscore -ii SequencescorerOne.txt "\
      "SequencescorerTwo.txt -k 4"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt7_k4.out"
  end
  
Name "gt sequencescorer Fscore TransProt11"
Keywords "gt_sequencescorer scorer fscore TransProt11"
Test do
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}SequencescorerOne.txt"
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}SequencescorerTwo.txt"  
  run "#{$bin}gt scorer -fscore -ii SequencescorerOne.txt "\
     "SequencescorerTwo.txt -k 5"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt11_k5.out"
end 


Name "gt sequencescorer Qgramdistance lowerbound unitcost"
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
  run "#{$bin}gt scorer -qgram -ii SequencescorerOne.txt "\
      "SequencescorerTwo.txt -q 4"
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

=begin
Name "gt sequencescorer Qgramdistance lowerbound blosom62"
Keywords "gt_sequencescorer scorer qgram lowerbound"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerOne.txt"
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}SequencescorerTwo.txt"
  run "#{$bin}gt scorer -qgram -ii SequencescorerOne.txt "\
      "SequencescorerTwo.txt -q 4"
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
  run "#{$bin}gt scorer -edist -ii SequencescorerOne.txt SequencescorerTwo.txt "\
      "-smatrix #{$testdata}BLOSUM62.out -indelscore -1"
  outfile = File.new("#{last_stdout}", "r")
  edistarray = outfile.readlines
  outfile.close
  for i in 0...edistarray.length
    if(edistarray[i][">"])
        next
    end
    tmp = (edistarray[i].split)
    edistarray.delete(edistarray[i])
    edistarray.insert(i,tmp[tmp.length-1])
  end
  for i in 0...edistarray.length
    if(("#{edistarray[i]}".to_i) < ("#{qgramarray[i]}".to_i)/8)
      #Rethink! in what relation are the score an the qgramdistance? 
    end
  end
end 
=end

Name "gt sequencescorer MaxMatches lowerbound unitcost"
Keywords "gt_sequencescorer scorer maxmatches lowerbound"
Test do
  edist = File.new("calcEdist.txt", "w")
  infile = File.new("#{$testdata}SequencescorerOne.txt", "r")
  inarray = infile.readlines
  infile.close
  for j in 0...inarray.length
    if(inarray[j][">"])
      next
    end
    break
  end
  for i in j+1...inarray.length
    if(inarray[i][">"])
      next
    end
    run "#{$bin}gt dev paircmp -e -ss #{inarray[j].chomp} #{inarray[i].chomp}"
    edist << last_stdout << "\n"
  end
  edist.close
  firstseq = File.new("FirstSeq.txt", "w")
  firstseq << inarray[0]
  firstseq << inarray[1] << "\n"
  firstseq.close
  run "#{$bin}gt suffixerator -db FirstSeq.txt -suf -tis -smap TransProt7 "\
      "-des no -md5 no -sds no"
  run "#{$bin}gt scorer -maxmatches -ii FirstSeq.txt "\
      "-seq #{$testdata}SequencescorerOne.txt"
  maxmatches = File.new("#{last_stdout}", "r")
  mmarray = maxmatches.readlines
  maxmatches.close
  for i in 0...mmarray.length
    if(mmarray[i][">"])
        next
    end
    if(mmarray[i]["#"])
        next
    end
    tmp = (mmarray[i].split)
    mmarray.delete(mmarray[i])
    mmarray.insert(i,tmp[tmp.length-1])
  end
  edist = File.new("calcEdist.txt", "r")
  edistarray = edist.readlines
  edist.close
  for i in 0...edistarray.length
    tmp = File.new("#{edistarray[i].chomp}", "r")
    tmparray = tmp.read.split()
    if(!(("#{mmarray[i]}".to_i) <= "#{tmparray[0]}".to_i))
      exit 1
    end
  end
end
