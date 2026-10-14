/**
 * @author Gabriel Henrique Silva
 * @file DaeSolution.cpp
 * @brief Solution-vector view implementation.
 */

#include "DaeSolution.h"
#include "Circuit.h"
#include "Node.h"

namespace e_XSpice {

int DaeSolution::rowOf(const Node* node) const {
    if (!node || node->isGround()) return -1;
    const int id = node->getId();
    if (id < 0 || static_cast<size_t>(id) >= node_map_.size()) return -1;
    return node_map_[id];
}

} // namespace e_XSpice
