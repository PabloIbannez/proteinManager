shared_proteinManager: dataTypes structure model chain residue atom
	g++-8 -shared -o libproteinManager.so dataTypes.o structure.o model.o chain.o residue.o atom.o
	rm *.o

dataTypes:
	g++-8 -fPIC -O3 -std=c++17 -I ../include -I ../include/proteinManager -c ../src/dataTypes.cpp -o dataTypes.o

structure:
	g++-8 -fPIC -O3 -std=c++17 -I ../include -I ../include/proteinManager -c ../src/structure.cpp -o structure.o

model:
	g++-8 -fPIC -O3 -std=c++17 -I ../include -I ../include/proteinManager -c ../src/model.cpp -o model.o

chain:
	g++-8 -fPIC -O3 -std=c++17 -I ../include -I ../include/proteinManager -c ../src/chain.cpp -o chain.o

residue:
	g++-8 -fPIC -O3 -std=c++17 -I ../include -I ../include/proteinManager -c ../src/residue.cpp -o residue.o

atom:
	g++-8 -fPIC -O3 -std=c++17 -I ../include -I ../include/proteinManager -c ../src/atom.cpp -o atom.o
