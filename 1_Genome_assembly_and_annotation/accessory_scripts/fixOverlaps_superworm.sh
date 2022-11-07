samtools faidx superworm_masked.fasta

cp test_removeContained_20.gff test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 549661)
loc3=$(echo 551124)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2\t551124/$var3\t551124/g" test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 565932)
loc3=$(echo 566281)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/565464\t$loc3/565464\t$var3/g" test_removeContained_21.gff
sed -i "s/565932\t$loc3/565932\t$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 2060878)
loc3=$(echo 2061217)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 17027745)
loc3=$(echo 17028467)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 18496966)
loc3=$(echo 18499644)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 21305364)
loc3=$(echo 21305714)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 30205658)
loc3=$(echo 30206065)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff
grep -v -w '30206117' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff
sed -i 's/30206149/30205946/g' test_removeContained_21.gff

loc1=$(echo scaffold_1)
loc2=$(echo 51121041)
loc3=$(echo 51121609)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 25774)
loc3=$(echo 26370)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2\t26370/$var3\t26370/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 7571598)
loc3=$(echo 7571824)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 12872891)
loc3=$(echo 12874501)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 12874353)
loc3=$(echo 12876134)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 44411981)
loc3=$(echo 44412251)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff
grep -v -w '44411763' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff
sed -i "s/44411702/44412154/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 48246407)
loc3=$(echo 48247879)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 59784093)
loc3=$(echo 59784416)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 61211706)
loc3=$(echo 61212556)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 61322795)
loc3=$(echo 61323174)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 62307131)
loc3=$(echo 62307538)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 65067012)
loc3=$(echo 65067222)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff
grep -w -v '65067270' test_removeContained_21.gff | grep -w -v '65071957' | grep -w -v '65072187' | grep -w -v '65077261' | grep -w -v '65077567' | grep -w -v '65077805' | grep -w -v '65078022' | grep -w -v '65078196' | grep -w -v '65085368' | grep -w -v '65085685' > temp.gff
sed -i 's/65086095/65067030/g' temp.gff
mv temp.gff test_removeContained_21.gff
grep -w 'Zmor_g8122' test_removeContained_20.gff | grep -v -w '65061040' | grep -v -w '65061372' | grep -v -w '65061547' | grep -v -w '65067222' | sed 's/65060500/65067270/g' | sed 's/gene7997/gene28594/g' | sed 's/Zmor_g8122/Zmor_g29044/g' | sed 's/mRNA8494/mRNA30463/g' >> test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 66198291)
loc3=$(echo 66198714)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_2)
loc2=$(echo 72828521)
loc3=$(echo 72828800)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_3)
loc2=$(echo 3135)
loc3=$(echo 3338)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2\t5828/$var3\t5828/g" test_removeContained_21.gff
sed -i "s/$loc2\t3338/$var3\t3338/g" test_removeContained_21.gff

loc1=$(echo scaffold_3)
loc2=$(echo 1615005)
loc3=$(echo 1615451)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_3)
loc2=$(echo 1714968)
loc3=$(echo 1715164)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_3)
loc2=$(echo 2163213)
loc3=$(echo 2163471)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_3)
loc2=$(echo 23676029)
loc3=$(echo 23676284)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_3)
loc2=$(echo 24718473)
loc3=$(echo 24718559)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
grep -v -w 'mRNA11280' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff

loc1=$(echo scaffold_3)
loc2=$(echo 38438980)
loc3=$(echo 38440800)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_4)
loc2=$(echo 4703981)
loc3=$(echo 4704403)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_4)
loc2=$(echo 7804294)
loc3=$(echo 7804491)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_4)
loc2=$(echo 26814430)
loc3=$(echo 26815130)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_4)
loc2=$(echo 46798900)
loc3=$(echo 46800192)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_4)
loc2=$(echo 56723440)
loc3=$(echo 56724198)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_5)
loc2=$(echo 1436848)
loc3=$(echo 1437243)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_5)
loc2=$(echo 4112802)
loc3=$(echo 4113276)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_5)
loc2=$(echo 5926655)
loc3=$(echo 5927216)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_5)
loc2=$(echo 14944590)
loc3=$(echo 14944883)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_5)
loc2=$(echo 20612606)
loc3=$(echo 20613898)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_5)
loc2=$(echo 30010115)
loc3=$(echo 30010470)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 2692397)
loc3=$(echo 2692777)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 3382677)
loc3=$(echo 3384017)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 5142082)
loc3=$(echo 5142672)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 9797860)
loc3=$(echo 9798340)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff
grep -w -v '9790190' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff
sed -i "s/9790157/9798031/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 19063411)
loc3=$(echo 19063581)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 23774618)
loc3=$(echo 23774869)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 25640557)
loc3=$(echo 25641777)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 33412935)
loc3=$(echo 33413219)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 40103372)
loc3=$(echo 40103677)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 41031416)
loc3=$(echo 41032030)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 45376173)
loc3=$(echo 45377473)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 545576)
loc3=$(echo 546453)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 4835803)
loc3=$(echo 4836080)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2\t4837270/$var3\t4837270/g" test_removeContained_21.gff
sed -i "s/$loc2\t4836080/$var3\t4836080/g" test_removeContained_21.gff
grep -v -w 'mRNA24077' test_removeContained_21.gff > temp.gff
sed -i 's/4837270/4828367/g' temp.gff
grep -w 'Zmor_g22980' test_removeContained_21.gff | grep -v 'mRNA24076' | sed 's/4826475/4835937/g' | sed 's/gene22636/gene28595/g' | sed 's/Zmor_g22980/Zmor_g29045/g' | sed 's/mRNA24077/mRNA30464/g' >> temp.gff
mv temp.gff test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 24036997)
loc3=$(echo 24037166)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 26145336)
loc3=$(echo 26146157)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 32059610)
loc3=$(echo 32060109)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 32578673)
loc3=$(echo 32578840)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 32633022)
loc3=$(echo 32633250)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 32956634)
loc3=$(echo 32956830)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 33048214)
loc3=$(echo 33048406)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 35516357)
loc3=$(echo 35516701)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 35559157)
loc3=$(echo 35559960)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 37589212)
loc3=$(echo 37589405)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_8)
loc2=$(echo 816595)
loc3=$(echo 817066)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_8)
loc2=$(echo 1050404)
loc3=$(echo 1050968)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_8)
loc2=$(echo 4639632)
loc3=$(echo 4640081)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_8)
loc2=$(echo 8284471)
loc3=$(echo 8285056)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_8)
loc2=$(echo 9837699)
loc3=$(echo 9838007)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_8)
loc2=$(echo 17436755)
loc3=$(echo 17437034)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_8)
loc2=$(echo 24011442)
loc3=$(echo 24012606)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 353592)
loc3=$(echo 354413)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/353592\t$loc3/353592\t$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 2673074)
loc3=$(echo 2673541)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 2937698)
loc3=$(echo 2939533)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 6650029)
loc3=$(echo 6650286)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 7933841)
loc3=$(echo 7934305)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 14168376)
loc3=$(echo 14168948)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 14360344)
loc3=$(echo 14360858)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 19662150)
loc3=$(echo 19662971)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_10)
loc2=$(echo 5364516)
loc3=$(echo 5365016)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_10)
loc2=$(echo 12533889)
loc3=$(echo 12534747)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff
grep -v -w '12533836' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff
sed -i "s/12533740/12534046/g" test_removeContained_21.gff

loc1=$(echo scaffold_10)
loc2=$(echo 16320989)
loc3=$(echo 16321576)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff
grep -w -v '16321658' test_removeContained_21.gff | grep -w -v '16324582' | grep -w -v '16324879' > temp.gff
mv temp.gff test_removeContained_21.gff
sed -i 's/16325656/16320580/g' test_removeContained_21.gff

loc1=$(echo scaffold_4)
loc2=$(echo 2139477)
loc3=$(echo 2139771)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff
grep -w -v '2138085' test_removeContained_21.gff > temp.gff
grep -w -v '2138712' temp.gff > test_removeContained_21.gff
grep -w -v '2138981' test_removeContained_21.gff > temp.gff
grep -w -v '2139214' temp.gff > test_removeContained_21.gff
grep -w -v '2139427' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff
sed -i 's/2137980/2139744/g' test_removeContained_21.gff

loc1=$(echo scaffold_6)
loc2=$(echo 45576526)
loc3=$(echo 45578529)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 1958683)
loc3=$(echo 1958892)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
grep -w -v 'mRNA23811' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff

loc1=$(echo scaffold_7)
loc2=$(echo 29612183)
loc3=$(echo 29612619)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff
grep -w -v '29596400' test_removeContained_21.gff > temp.gff
grep -w -v '29607487' temp.gff > test_removeContained_21.gff
sed -i 's/29596247/29612414/g' test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 1702063)
loc3=$(echo 1702756)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/1702063\t$loc3/1702063\t$var3/g" test_removeContained_21.gff
grep -w -v '1702919' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff
sed -i 's/1703003/1702539/g' test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 8466052)
loc3=$(echo 8466663)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
grep -v -w 'Zmor_g28147' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 9155104)
loc3=$(echo 9155596)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f1,1 -d'N' | wc -c)-2))")
var2=$(echo $loc2)
var3=$(echo "$(($var2+$var))")
grep "$loc3" test_removeContained_21.gff
sed -i "s/$loc3/$var3/g" test_removeContained_21.gff
grep -w -v '9156658' test_removeContained_21.gff > temp.gff
grep -w -v '9156890' temp.gff > test_removeContained_21.gff
sed -i 's/9157198/9155493/g' test_removeContained_21.gff

loc1=$(echo scaffold_9)
loc2=$(echo 16254546)
loc3=$(echo 16254884)
samtools faidx superworm_masked.fasta "$loc1:$loc2-$loc3" | grep -v '>'
grep -w -v '16254478' test_removeContained_21.gff > temp.txt
mv temp.txt test_removeContained_21.gff
var=$(echo "$(($(samtools faidx superworm_masked.fasta $loc1:$loc2-$loc3 | grep -v '>' | tr --delete '\n' | cut -f101,101 -d'N' | wc -c)-2))")
var2=$(echo $loc3)
var3=$(echo "$(($var2-$var))")
grep "$loc2" test_removeContained_21.gff
sed -i "s/$loc2/$var3/g" test_removeContained_21.gff
sed -i "s/16254237/16254656/g" test_removeContained_21.gff

grep -w -v 'mRNA16951' test_removeContained_21.gff > temp.gff
sed -i 's/56724198/56721021/g' temp.gff
grep 'ID=gene15980' test_removeContained_21.gff | sed 's/gene15980/gene28596/g' | sed 's/Zmor_g16211/Zmor_g29046/g' | sed 's/56719696/56723544/g' >> temp.gff
grep -w 'Zmor_g16211' test_removeContained_21.gff | grep 'mRNA16951' | sed 's/gene15980/gene28596/g' | sed 's/Zmor_g16211/Zmor_g29046/g' | sed 's/mRNA16951/mRNA30468/g' >> temp.gff
mv temp.gff test_removeContained_21.gff

grep -w -v 'mRNA17424' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff
sed -i 's/7383\t39176/35943\t39176/g' test_removeContained_21.gff
sed -i 's/7383 39176/35943 39176/g' test_removeContained_21.gff

grep -w '17844075' test_removeContained_21.gff
grep -v -w 'Zmor_g842' test_removeContained_21.gff > temp.gff
grep -w '41524293' test_removeContained_21.gff
grep -v -w 'Zmor_g2132' temp.gff > test_removeContained_21.gff
grep -w '54097519' test_removeContained_21.gff
grep -v -w 'Zmor_g2628' test_removeContained_21.gff > temp.gff
grep -w '13627981' test_removeContained_21.gff
grep -v -w 'Zmor_g5251' temp.gff > test_removeContained_21.gff
grep -w '13759997' test_removeContained_21.gff
grep -v -w 'Zmor_g5256' test_removeContained_21.gff > temp.gff
grep -w '26925184' test_removeContained_21.gff
grep -v -w 'Zmor_g5846' temp.gff > test_removeContained_21.gff
grep -w '58610037' test_removeContained_21.gff
grep -v -w 'Zmor_g7722' test_removeContained_21.gff > temp.gff
grep -w '58616887' test_removeContained_21.gff
grep -v -w 'Zmor_g7723' temp.gff > test_removeContained_21.gff
grep -w '58622051' test_removeContained_21.gff
grep -v -w 'Zmor_g7724' test_removeContained_21.gff > temp.gff
grep -w '68711465' test_removeContained_21.gff
grep -v -w 'Zmor_g8332' temp.gff > test_removeContained_21.gff
grep -w '9302831' test_removeContained_21.gff
grep -v -w 'Zmor_g9864' test_removeContained_21.gff > temp.gff
grep -w '13206735' test_removeContained_21.gff
grep -v -w 'Zmor_g10218' temp.gff > test_removeContained_21.gff
grep -w '28759348' test_removeContained_21.gff
grep -v -w 'Zmor_g10968' test_removeContained_21.gff > temp.gff
grep -w '35645029' test_removeContained_21.gff
grep -v -w 'Zmor_g11331' temp.gff > test_removeContained_21.gff
grep -w '35760513' test_removeContained_21.gff
grep -v -w 'Zmor_g11340' test_removeContained_21.gff > temp.gff
grep -w '2020521' test_removeContained_21.gff
grep -v -w 'Zmor_g12849' temp.gff > test_removeContained_21.gff
grep -w '11204279' test_removeContained_21.gff
grep -v -w 'Zmor_g13453' test_removeContained_21.gff > temp.gff
grep -w '11773981' test_removeContained_21.gff
grep -v -w 'Zmor_g13484' temp.gff > test_removeContained_21.gff
grep -w '29729950' test_removeContained_21.gff
grep -v -w 'Zmor_g14414' test_removeContained_21.gff > temp.gff
grep -w '29810037' test_removeContained_21.gff
grep -v -w 'Zmor_g14416' temp.gff > test_removeContained_21.gff
grep -w '32736562' test_removeContained_21.gff
grep -v -w 'Zmor_g14611' test_removeContained_21.gff > temp.gff
grep -w '38366877' test_removeContained_21.gff
grep -v -w 'Zmor_g14960' temp.gff > test_removeContained_21.gff
grep -w '40713396' test_removeContained_21.gff
grep -v -w 'Zmor_g15080' test_removeContained_21.gff > temp.gff
grep -w '40766113' test_removeContained_21.gff
grep -v -w 'Zmor_g15082' temp.gff > test_removeContained_21.gff
grep -w '42272809' test_removeContained_21.gff
grep -v -w 'Zmor_g15184' test_removeContained_21.gff > temp.gff
grep -w '766411' test_removeContained_21.gff
grep -v -w 'Zmor_g16825' temp.gff > test_removeContained_21.gff
grep -w '2070911' test_removeContained_21.gff
grep -v -w 'Zmor_g16991' test_removeContained_21.gff > temp.gff
grep -w '2097026' test_removeContained_21.gff
grep -v -w 'Zmor_g16992' temp.gff > test_removeContained_21.gff
grep -w '2258488' test_removeContained_21.gff
grep -v -w 'Zmor_g17000' test_removeContained_21.gff > temp.gff
grep -w '2258992' test_removeContained_21.gff
grep -v -w 'Zmor_g17001' temp.gff > test_removeContained_21.gff
grep -w '3942527' test_removeContained_21.gff
grep -v -w 'Zmor_g17090' test_removeContained_21.gff > temp.gff
grep -w '4052847' test_removeContained_21.gff
grep -v -w 'Zmor_g17091' temp.gff > test_removeContained_21.gff
grep -w '21941597' test_removeContained_21.gff
grep -v -w 'Zmor_g18422' test_removeContained_21.gff > temp.gff
grep -w '8358864' test_removeContained_21.gff
grep -v -w 'Zmor_g20108' temp.gff > test_removeContained_21.gff
grep -w '18686366' test_removeContained_21.gff
grep -v -w 'Zmor_g20870' test_removeContained_21.gff > temp.gff
grep -w '33972381' test_removeContained_21.gff
grep -v -w 'Zmor_g21705' temp.gff > test_removeContained_21.gff
grep -w '35767696' test_removeContained_21.gff
grep -v -w 'Zmor_g21797' test_removeContained_21.gff > temp.gff
grep -w '38950391' test_removeContained_21.gff
grep -v -w 'Zmor_g21955' temp.gff > test_removeContained_21.gff
grep -w '3874375' test_removeContained_21.gff
grep -v -w 'Zmor_g22890' test_removeContained_21.gff > temp.gff
grep -w '10413231' test_removeContained_21.gff
grep -v -w 'Zmor_g23391' temp.gff > test_removeContained_21.gff
grep -w '15890988' test_removeContained_21.gff
grep -v -w 'Zmor_g23638' test_removeContained_21.gff > temp.gff
grep -w '18045559' test_removeContained_21.gff
grep -v -w 'Zmor_g23770' temp.gff > test_removeContained_21.gff
grep -w '19331622' test_removeContained_21.gff
grep -v -w 'Zmor_g23828' test_removeContained_21.gff > temp.gff
grep -w '22319602' test_removeContained_21.gff
grep -v -w 'Zmor_g24039' temp.gff > test_removeContained_21.gff
grep -w '26578650' test_removeContained_21.gff
grep -v -w 'Zmor_g24279' test_removeContained_21.gff > temp.gff
grep -w '26596733' test_removeContained_21.gff
grep -v -w 'Zmor_g24280' temp.gff > test_removeContained_21.gff
grep -w '5307168' test_removeContained_21.gff
grep -v -w 'Zmor_g25651' test_removeContained_21.gff > temp.gff
grep -w '21493672' test_removeContained_21.gff
grep -v -w 'Zmor_g26677' temp.gff > test_removeContained_21.gff
grep -w '23892260' test_removeContained_21.gff
grep -v -w 'Zmor_g26822' test_removeContained_21.gff > temp.gff
grep -w '29964345' test_removeContained_21.gff
grep -v -w 'Zmor_g27206' temp.gff > test_removeContained_21.gff
grep -w '27416' test_removeContained_21.gff
grep -v -w 'Zmor_g16683' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff

sed -i 's/2085\t5775/1993\t5775/' test_removeContained_21.gff
sed -i 's/2085\t3544/1993\t3544/' test_removeContained_21.gff
sed -i 's/2085\t2247/1993\t2247/' test_removeContained_21.gff

grep -w -v '3910547' test_removeContained_21.gff > temp.gff
mv temp.gff test_removeContained_21.gff
sed -i 's/3865360/3865402/' test_removeContained_21.gff
sed -i 's/3910699/3865402/' test_removeContained_21.gff
