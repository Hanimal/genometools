Name "gt_sequencescorer_check_gramdistance"
Keywords "gt_sequencescorer qgram distance"
Test do
  run "#{$bin}gt scorer -qgram -s #{$testdata}sequencescorer_one.txt #{$testdata}sequencescorer_two.txt -q 4"
  run "diff -B #{last_stdout} #{$testdata}sequencescorer_compare_q=4.out"
  run "#{$bin}gt scorer -qgram -s #{$testdata}sequencescorer_one.txt #{$testdata}sequencescorer_two.txt -q 3"
  run "diff -B #{last_stdout} #{$testdata}sequencescorer_compare_q=3.out"
end 
