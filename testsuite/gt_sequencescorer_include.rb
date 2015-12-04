class FastaIterator
  def initialize(filename)
    @file = get_file(filename)
  end

  def each
    seq_list = Array.new()
    header = nil
    @file.each_line do |line|
      if not line.match(/^(>|\s*$|\s*;)/)
        seq_list.push(line.chomp)
      elsif m = line.match(/^>(.*)/)
        if header != nil
          yield header, seq_list.join("").gsub(/\s/,"")
          seq_list.clear
        end
        header = m[1].chomp
      end
    end
    if header != nil
      yield header, seq_list.join("").gsub(/\s/,"")
    end
    @file.rewind
  end

  private

  def get_file(filename)
    file = nil
    if File.exist?(filename)
      begin
        file = File.new(filename,"r")
      rescue => err
        STDERR.puts "Could not open file #{filename}: #{err}"
        exit 1
      end
    else
      STDERR.puts "File #{filename} does not exist!"
      exit 1
    end
    return file
  end
end

Name "gt sequencescorer Alphabet"
Keywords "gt_sequencescorer scorer alphabet"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}sw100K1.fsa"
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}sw100K2.fsa"
  run "#{$bin}gt scorer -qgram -ii sw100K1.fsa "\
      " sw100K2.fas -q 4", :retval => 1
end

Name "gt sequencescorer Fileoption"
Keywords "gt_sequencescorer scorer file"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}sw100K1.fsa"
  run "#{$bin}gt scorer -qgram -ii sw100K1.fsa -q 4 -file", :retval => 1
  run "#{$bin}gt scorer -fscore -ii sw100K1.fsa -k 4 -file", :retval => 1
end

Name "gt sequencescorer Fileformat"
Keywords "gt_sequencescorer scorer Fileformat"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}sw100K1.fsa"
  run "#{$bin}gt scorer -qgram -ii sw100K1.fsa "\
      "sw100K2.fas -q 4", :retval => 1
  run "#{$bin}gt scorer -qgram -ii -q 4", :retval => 1
  run "#{$bin}gt scorer -qgram -ii sw100K1.fsa "\
      "sw100K1.fas sw100K1.fsa -q 4", :retval => 1
end

Name "gt sequencescorer Qgramdistance"
Keywords "gt_sequencescorer scorer qgram"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}sw100K1.fsa"
  run "#{$bin}gt scorer -qgram -ii sw100K1.fsa -q 4", :retval => 0
end

Name "gt sequencescorer Qgramdistance TransProt7"
Keywords "gt_sequencescorer scorer qgram"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}sw100K1.fsa"
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}sw100K2.fsa"
  run "#{$bin}gt scorer -qgram -ii sw100K1.fsa "\
      "sw100K2.fsa -q 4 -distance"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt7_q4.out"
end

Name "gt sequencescorer Qgramdistance TransProt11"
Keywords "gt_sequencescorer scorer qgram TransProt11"
Test do
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}sw100K1.fsa"
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}sw100K2.fsa"
  run "#{$bin}gt scorer -qgram -ii sw100K1.fsa "\
      "sw100K2.fsa -q 5 -distance"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt11_q5.out"
end

Name "gt sequencescorer Fscore"
Keywords "gt_sequencescorer scorer fscore"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}sw100K1.fsa"
  run_test "#{$bin}gt scorer -fscore -ii sw100K1.fsa -k 6", :retval => 0
end

Name "gt sequencescorer Fscore TransProt7"
Keywords "gt_sequencescorer scorer fscore"
Test do
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}sw100K1.fsa"
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}sw100K2.fsa"
  run "#{$bin}gt scorer -fscore -ii sw100K1.fsa "\
      "sw100K2.fsa -k 4"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt7_k4.out"
  end
  
Name "gt sequencescorer Fscore TransProt11"
Keywords "gt_sequencescorer scorer fscore"
Test do
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}sw100K1.fsa"
  run "#{$bin}gt encseq encode -smap TransProt11 #{$testdata}sw100K2.fsa"  
  run "#{$bin}gt scorer -fscore -ii sw100K1.fsa "\
     "sw100K2.fsa -k 5"
  run "diff -B #{last_stdout} #{$testdata}sequencescorerCompare_TransProt11_k5.out"
end 

Name "gt sequencescorer Maxmatches TransProt7"
Keywords "gt_sequencescorer scorer maxmatches"
Test do
  outfile = File.new("maxmatchesout","w")
  idx = 0
  `#{$bin}gt encseq encode -smap TransProt7 #{$testdata}sw100K1.fsa`
  sequencefile = "#{$testdata}sw100K1.fsa".split("/")
  sequencefile = sequencefile[sequencefile.length-1]
  data = File.new("data.fas", "w")
  data << `#{$bin}gt encseq decode #{sequencefile}`
  data.close
    
  fi = FastaIterator.new("data.fas")
  fi.each do |header,sequence|
    tmp = File.new("tmp.fas", "w")
    tmp << ">" << header << "\n"
    tmp << sequence << "\n"
    tmp.close
    `#{$bin}gt suffixerator -db tmp.fas -suf -tis -des no -md5 no -sds no -smap TransProt7`
    maxmatch = `#{$bin}gt scorer -maxmatches -ii tmp.fas -seq data.fas -distance`
    maxmatch = maxmatch.split
    for i in 0...maxmatch.length
      if(maxmatch[i] == "Maxmatchesdistance" or maxmatch[i] == "in" or \
         maxmatch[i] == "sequence" or maxmatch[i] == "is")
        maxmatch.delete(maxmatch[i])
      end
    end
    i = 1
    while(i < maxmatch.length)
      if(idx > maxmatch[i-1].to_i)
        outfile << "Maxmatches from #{maxmatch[i-1]} in #{idx} is #{maxmatch[i]}\n"
      end
      i = i + 2
    end
    idx = idx + 1
    `rm tmp.*`    
  end
  `rm data.fas`
  outfile.close

  run "diff -B maxmatchesout #{$testdata}sequencescorerCompare_TransProt7_maxmatches"
end 

Name "gt sequencescorer Qgramdistance lowerbound"
Keywords "gt_sequencescorer scorer qgram lowerbound"
Test do
  outfile = File.new("edist.txt", "w")
  fi1 = FastaIterator.new("#{$testdata}sw100K1.fsa")
  fi2 = FastaIterator.new("#{$testdata}sw100K2.fsa")
  fi1.each do |header1,sequence1|
    fi2.each do |header2,sequence2|
        run "#{$bin}gt dev paircmp -e -ss #{sequence1} #{sequence2}"
        outfile << last_stdout << "\n"
    end
  end
  outfile.close
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}sw100K1.fsa"
  run "#{$bin}gt encseq encode -smap TransProt7 #{$testdata}sw100K2.fsa"
  run "#{$bin}gt scorer -qgram -ii sw100K1.fsa sw100K2.fsa -q 4 -distance"
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
  outfile = File.new("edist.txt", "r")
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

Name "gt sequencescorer MaxMatches lowerbound"
Keywords "gt_sequencescorer scorer maxmatches lowerbound"
Test do
  tmp = ""
  edist = File.new("calcEdist.txt", "w")
  `#{$bin}gt encseq encode -smap TransProt7 #{$testdata}sw100K1.fsa`
  data = File.new("data.fas", "w")
  data << `#{$bin}gt encseq decode sw100K1.fsa`
  data.close
  fi = FastaIterator.new("data.fas")
  fi.each do |header,sequence|
    tmp = sequence
    break
  end
  fi.each do |header,sequence|
    run "#{$bin}gt dev paircmp -e -ss #{tmp} #{sequence}"
    edist << last_stdout << "\n"
  end
  edist.close
  firstseq = File.new("FirstSeq.txt", "w")
  firstseq << ">FirstSeq \n" << tmp << "\n" 
  firstseq.close
  run "#{$bin}gt suffixerator -db FirstSeq.txt -suf -tis -smap TransProt7 "\
      "-des no -md5 no -sds no"
  run "#{$bin}gt scorer -maxmatches -ii FirstSeq.txt "\
      "-seq data.fas -distance"
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
