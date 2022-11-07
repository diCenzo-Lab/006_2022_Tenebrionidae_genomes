samtools faidx mealworm_masked.fasta

cp test_removeContained_20.gff test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 1075285)
loc3=$(echo 1075695)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 2527362)
loc3=$(echo 2527850)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 21090372)
loc3=$(echo 21091160)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 22831473)
loc3=$(echo 22832009)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 32195024)
loc3=$(echo 32195563)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 4338612)
loc3=$(echo 4339469)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 10042755)
loc3=$(echo 10043444)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_3)
loc2=$(echo 4949919)
loc3=$(echo 4951499)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_4)
loc2=$(echo 9649424)
loc3=$(echo 9649661)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
grep -v -w '9649147' test_removeContained_21.gff > temp.gff
grep -v -w '9649363' temp.gff > test_removeContained_21.gff
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff
sed -i "s/9648963/9649535/g" test_removeContained_21.gff

loc1=$(echo scaffold_4)
loc2=$(echo 15029246)
loc3=$(echo 15029826)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_4)
loc2=$(echo 15913041)
loc3=$(echo 15913314)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff
grep -v -w '15921661' test_removeContained_21.gff > temp.gff
grep -v -w '15921699' temp.gff > test_removeContained_21.gff
sed -i "s/15921714/15913314/g" test_removeContained_21.gff

loc1=$(echo scaffold_4)
loc2=$(echo 18156568)
loc3=$(echo 18157413)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_5)
loc2=$(echo 20738981)
loc3=$(echo 20740132)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 2049723)
loc3=$(echo 2050436)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 19895339)
loc3=$(echo 19898407)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 20145525)
loc3=$(echo 20146061)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 6728983)
loc3=$(echo 6729318)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
grep -v -w 'Tmol_g17084' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 10865579)
loc3=$(echo 10866022)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 16862693)
loc3=$(echo 16864033)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 3261188)
loc3=$(echo 3261546)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 4329326)
loc3=$(echo 4329578)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
grep -v -w 'mRNA20810' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 11418365)
loc3=$(echo 11418985)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_11)
loc2=$(echo 4936676)
loc3=$(echo 4939096)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_12)
loc2=$(echo 339776)
loc3=$(echo 340192)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/339776\t$loc3/339776\t$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_12)
loc2=$(echo 4249035)
loc3=$(echo 4249937)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
grep -v -w 'mRNA4941' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff

loc1=$(echo scaffold_12)
loc2=$(echo 5057179)
loc3=$(echo 5057595)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_12)
loc2=$(echo 5876962)
loc3=$(echo 5877879)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_13)
loc2=$(echo 4836150)
loc3=$(echo 4836984)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_13)
loc2=$(echo 6348055)
loc3=$(echo 6348391)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
grep -v -w 'mRNA5898' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff

loc1=$(echo scaffold_13)
loc2=$(echo 6366086)
loc3=$(echo 6366488)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
grep -v -w 'mRNA5900' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff

loc1=$(echo scaffold_14)
loc2=$(echo 676843)
loc3=$(echo 679077)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_14)
loc2=$(echo 1631762)
loc3=$(echo 1635889)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff
grep -v -w '1636379' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff
sed -i "s/1636924/1635786/g" test_removeContained_21.gff

loc1=$(echo scaffold_14)
loc2=$(echo 3135302)
loc3=$(echo 3135688)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$var3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2\t3135688/$var3\t3135688/g" test_removeContained_21.gff

loc1=$(echo scaffold_14)
loc2=$(echo 4177878)
loc3=$(echo 4182143)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_15)
loc2=$(echo 1432477)
loc3=$(echo 1432917)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_15)
loc2=$(echo 1330057)
loc3=$(echo 1330615)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff
grep -v -w '1330676' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff
sed -i "s/1330754/1330388/g" test_removeContained_21.gff

loc1=$(echo scaffold_15)
loc2=$(echo 3366700)
loc3=$(echo 3367040)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_16)
loc2=$(echo 2584692)
loc3=$(echo 2586608)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_16)
loc2=$(echo 2619906)
loc3=$(echo 2620934)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 1304492)
loc3=$(echo 1305052)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_5)
loc2=$(echo 20401694)
loc3=$(echo 20402143)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
grep -v -w 'mRNA15479' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 7439719)
loc3=$(echo 7440021)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff
grep -v -w '7438214' test_removeContained_21.gff > temp.gff
grep -v -w '7438161' temp.gff > test_removeContained_21.gff
sed -i "s/7438072/7439901/g" test_removeContained_21.gff

loc1=$(echo scaffold_12)
loc2=$(echo 2139613)
loc3=$(echo 2140333)
samtools faidx mealworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx mealworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

grep -v -w 'Tmol_g8243' test_removeContained_21.gff > temp.gff
grep -v -w 'Tmol_g13029' temp.gff > test_removeContained_21.gff
grep -v -w 'mRNA8349' test_removeContained_21.gff > temp.gff
grep -v -w 'mRNA5016' temp.gff > test_removeContained_21.gff

grep -v -w 'mRNA1995' test_removeContained_21.gff > temp.gff
sed -i 's/22832102/22826972/g' temp.gff
grep -w 'Tmol_g1845' test_removeContained_21.gff | sed 's/gene1811/gene19870/g' | sed 's/Tmol_g1845/Tmol_g20219/g' | sed 's/mRNA1995/mRNA21455/g' | grep -v -w 'mRNA1994' | sed 's/22826523/22831596/g' >> temp.gff
mv temp.gff test_removeContained_21.gff

grep -v -w 'mRNA7500' test_removeContained_21.gff > temp.gff
sed -i 's/2586608/2583337/g' temp.gff
grep -w 'Tmol_g6999' test_removeContained_21.gff | sed 's/gene6870/gene19871/g' | sed 's/Tmol_g6999/Tmol_g20220/g' | sed 's/mRNA7500/mRNA21456/g' | grep -v -w 'mRNA7499' | sed 's/2582714/2584807/g' >> temp.gff
mv temp.gff test_removeContained_21.gff

grep -w '3720997' test_removeContained_21.gff
grep -v -w 'Tmol_g545' test_removeContained_21.gff > temp.gff
grep -w '3906036' test_removeContained_21.gff
grep -v -w 'Tmol_g596' temp.gff > test_removeContained_21.gff
grep -w '7112904' test_removeContained_21.gff
grep -v -w 'Tmol_g896' test_removeContained_21.gff > temp.gff
grep -w '16954187' test_removeContained_21.gff
grep -v -w 'Tmol_g1487' temp.gff > test_removeContained_21.gff
grep -w '24478057' test_removeContained_21.gff
grep -v -w 'Tmol_g1935' test_removeContained_21.gff > temp.gff
grep -w '24480896' test_removeContained_21.gff
grep -v -w 'Tmol_g1936' temp.gff > test_removeContained_21.gff
grep -w '73924' test_removeContained_21.gff
grep -v -w 'Tmol_g7633' test_removeContained_21.gff > temp.gff
grep -w '2811908' test_removeContained_21.gff
grep -v -w 'Tmol_g7721' temp.gff > test_removeContained_21.gff
grep -w '18624920' test_removeContained_21.gff
grep -v -w 'Tmol_g8561' test_removeContained_21.gff > temp.gff
grep -w '7529311' test_removeContained_21.gff
grep -v -w 'Tmol_g9788' temp.gff > test_removeContained_21.gff
grep -w '7529923' test_removeContained_21.gff
grep -v -w 'Tmol_g9789' test_removeContained_21.gff > temp.gff
grep -w '5382684' test_removeContained_21.gff
grep -v -w 'Tmol_g12248' temp.gff > test_removeContained_21.gff
grep -w '5741045' test_removeContained_21.gff
grep -v -w 'Tmol_g12315' test_removeContained_21.gff > temp.gff
grep -w '14151098' test_removeContained_21.gff
grep -v -w 'Tmol_g12894' temp.gff > test_removeContained_21.gff
grep -w '6378704' test_removeContained_21.gff
grep -v -w 'Tmol_g13984' test_removeContained_21.gff > temp.gff
grep -w '6384331' test_removeContained_21.gff
grep -v -w 'Tmol_g13985' temp.gff > test_removeContained_21.gff
grep -w '6413163' test_removeContained_21.gff
grep -v -w 'Tmol_g13986' test_removeContained_21.gff > temp.gff
grep -w '6886449' test_removeContained_21.gff
grep -v -w 'Tmol_g14031' temp.gff > test_removeContained_21.gff
grep -w '8762956' test_removeContained_21.gff
grep -v -w 'Tmol_g14133' test_removeContained_21.gff > temp.gff
grep -w '18319864' test_removeContained_21.gff
grep -v -w 'Tmol_g14513' temp.gff > test_removeContained_21.gff
grep -w '21487536' test_removeContained_21.gff
grep -v -w 'Tmol_g14631' test_removeContained_21.gff > temp.gff
grep -w '19267887' test_removeContained_21.gff
grep -v -w 'Tmol_g16184' temp.gff > test_removeContained_21.gff
grep -w '1702501' test_removeContained_21.gff
grep -v -w 'Tmol_g16538' test_removeContained_21.gff > temp.gff
grep -w '1703946' test_removeContained_21.gff
grep -v -w 'Tmol_g16539' temp.gff > test_removeContained_21.gff
grep -w '2025187' test_removeContained_21.gff
grep -v -w 'Tmol_g16598' test_removeContained_21.gff > temp.gff
grep -w '8642650' test_removeContained_21.gff
grep -v -w 'Tmol_g17202' temp.gff > test_removeContained_21.gff
grep -w '12133248' test_removeContained_21.gff
grep -v -w 'Tmol_g17411' test_removeContained_21.gff > temp.gff
grep -w '2144370' test_removeContained_21.gff
grep -v -w 'Tmol_g17937' temp.gff > test_removeContained_21.gff
grep -w '5429494' test_removeContained_21.gff
grep -v -w 'Tmol_g18345' test_removeContained_21.gff > temp.gff
grep -w '14919998' test_removeContained_21.gff
grep -v -w 'Tmol_g19200' temp.gff > test_removeContained_21.gff
grep -w '14986707' test_removeContained_21.gff
grep -v -w 'Tmol_g19207' test_removeContained_21.gff > temp.gff
grep -w '3414589' test_removeContained_21.gff
grep -v -w 'Tmol_g2625' temp.gff > test_removeContained_21.gff
grep -w '3415036' test_removeContained_21.gff
grep -v -w 'Tmol_g2626' test_removeContained_21.gff > temp.gff
grep -w '6631617' test_removeContained_21.gff
grep -v -w 'Tmol_g2835' temp.gff > test_removeContained_21.gff
grep -w '7754954' test_removeContained_21.gff
grep -v -w 'Tmol_g2897' test_removeContained_21.gff > temp.gff
grep -w '1151441' test_removeContained_21.gff
grep -v -w 'Tmol_g3281' temp.gff > test_removeContained_21.gff
grep -w '6158992' test_removeContained_21.gff
grep -v -w 'Tmol_g3682' test_removeContained_21.gff > temp.gff
grep -w '223366' test_removeContained_21.gff
grep -v -w 'Tmol_g7576' temp.gff > test_removeContained_21.gff
