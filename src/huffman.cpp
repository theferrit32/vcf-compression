#include <map>
#include <vector>

typedef uint8_t byte_t;

template <class T> class LinkedBinaryNode {
public:
    LinkedBinaryNode(T *data) {
        this->data = data;
    }

    LinkedBinaryNode<T> *left;
    LinkedBinaryNode<T> *right;
    T *data;
    size_t count;
};

class ByteHuffman {
public:
    /**
     *
     */
    ByteHuffman(std::vector<byte_t> data_vec) {
        for (auto iter = data_vec.begin(); iter != data_vec.end(); iter++) {
            byte_t d = *iter;

        }
    }
private:
    std::map<byte_t, LinkedBinaryNode<byte_t>*> data_frequency_map;


};


/**
 * Vector-based binary tree
 */
// template <class T> class UniqueBinaryTree {
// public:
//     UniqueBinaryTree() {
//         this->root = NULL;
//     }

//     /**
//      * T pointer must be heap-allocated or exist for the lifespan of the tree.
//      */
//     void union(LinkedBinaryNode) {
//         LinkedBinaryNode *node = new LinkedBinaryNode(obj);
//         if (root == NULL) {
//             root = node;
//             return;
//         }
//     }



// private:
//     LinkedBinaryNode<T> *root;
//     std::map<T, LinkedBinaryNode<T>> node_map;
// };

