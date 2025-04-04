# R
echo ""; echo "Checking for R..."; echo ""
R --version | head -n 1

# R PACKAGES
echo ""; echo "Testing for R packages..."; echo ""
UE="UE"; SE="SE"
for package in "Rcpp" "devtools" "R6" "ape" "phytools" "distributions3" "sf" "hms" "condMVNorm" "mcmcse" "rase" "spatstat"
do
    output=`echo $(R -e "\"${package}\" %in% rownames(installed.packages())") | tail -c 7 | head -c 2`
    if [ "$output" = "$UE" ]; then echo "$package installed"; else echo "$package not found"; fi
done
echo ""; echo "Full List:"; echo ""
echo $(R -e "cat(rownames(installed.packages()))") | cut -d ">" -f2- | cut -c 38- | rev | cut -c 4- | rev