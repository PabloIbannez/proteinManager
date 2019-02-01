shared_proteinManager: dataTypes structure model chain residue atom
	g++ -shared -o libproteinManager.so dataTypes.o structure.o model.o chain.o residue.o atom.o
	rm *.o

dataTypes:
	g++ -fPIC -O3 -std=c++11 -I ../include -c ../src/dataTypes.cpp -o dataTypes.o

structure:
	g++ -fPIC -O3 -std=c++11 -I ../include -c ../src/structure.cpp -o structure.o

model:
	g++ -fPIC -O3 -std=c++11 -I ../include -c ../src/model.cpp -o model.o

chain:
	g++ -fPIC -O3 -std=c++11 -I ../include -c ../src/chain.cpp -o chain.o

residue:
	g++ -fPIC -O3 -std=c++11 -I ../include -c ../src/residue.cpp -o residue.o

atom:
	g++ -fPIC -O3 -std=c++11 -I ../include -c ../src/atom.cpp -o atom.o
