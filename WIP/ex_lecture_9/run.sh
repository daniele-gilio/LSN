#----------------Remove Files-------------------
rm best*
rm mean*
#----------------Circumference Distribution----------------------
sed -i '3s/.*/150/' input.dat #Generation Number
sed -i '4s/.*/0/' input.dat #Mercy Number
sed -i '6s/.*/0/' input.dat #City Distribution
./tsp.exe
#----------------Square Distribution----------------------
sed -i '3s/.*/250/' input.dat #Generation Number
sed -i '4s/.*/8/' input.dat #Mercy Number
sed -i '6s/.*/1/' input.dat #City Distribution
./tsp.exe
