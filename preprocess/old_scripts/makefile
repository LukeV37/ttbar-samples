default: split_test_train.cpp
	g++ -o ./split split_test_train.cpp -std=c++11
	g++ -o ./edges create_edges.cpp -std=c++11

.PHONY: clean split edges

split:
	./split
edges:
	./edges

clean:
	rm split 
	rm edges
