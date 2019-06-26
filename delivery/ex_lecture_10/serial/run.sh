#----------------Remove Files-------------------
rm best*
#----------------Circumference Distribution----------------------
sed -i '2s/.*/10000/' input.dat #Generation Number
sed -i '4s/.*/0/' input.dat #City Distribution
./tsp.exe
#----------------Square Distribution----------------------
sed -i '2s/.*/10000/' input.dat #Generation Number
sed -i '4s/.*/1/' input.dat #City Distribution
./tsp.exe
