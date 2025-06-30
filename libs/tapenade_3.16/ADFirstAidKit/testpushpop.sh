cleantapstack() {
  ls tapStack* 2>/dev/null
  rm tapStack* 2>/dev/null
}

for i in {1..999}
 do
  echo "SEED:"$i
  ./testpushpop 0 "*" $i ; cleantapstack
  ./testpushpop 0 "3[* *]" $i ; cleantapstack
  ./testpushpop 0 "* 2[* 2[*] *]" $i ; cleantapstack
  ./testpushpop 0 "* * (* *) * *" $i ; cleantapstack
  ./testpushpop 0 "* * 2[* * * *] * *" $i ; cleantapstack
  ./testpushpop 0 "* * 2[* * (* * * *) *] * *" $i ; cleantapstack
  ./testpushpop 0 "* * 2[* * 3[*] *] * *" $i ; cleantapstack
  ./testpushpop 0 "* * 2[* (* 3[*] *) * (3[((*)*)*]) * * *] *" $i ; cleantapstack
  ./testpushpop 0 "* (*) * L(* *) *" $i ; cleantapstack
  ./testpushpop 0 "* L(* (* L((*) *) *) *) *" $i ; cleantapstack
  ./testpushpop 0 "* 3[* L(* *) *] *" $i ; cleantapstack
  ./testpushpop 0 "* (* 2[* 3[* L(* *) *] *] *) (*) *" $i ; cleantapstack
  ./testpushpop 0 "* 3[2[(*) *] (* *) *] *" $i ; cleantapstack
  ./testpushpop 0 "* 3[(* *) * 2[(*) *]] *" $i ; cleantapstack
  ./testpushpop 0 "* (((*) *) * (* (*))) *" $i ; cleantapstack
  ./testpushpop 0 "* (2[((*) *)] 3[* (* (*))]) *" $i ; cleantapstack
  ./testpushpop 0 "* 3[(2[*])]"  $i ; cleantapstack
  ./testpushpop 0 "* 3[(2[(*)])]"  $i ; cleantapstack
  ./testpushpop 0 "* 2[* L(* *)]"  $i ; cleantapstack
done
