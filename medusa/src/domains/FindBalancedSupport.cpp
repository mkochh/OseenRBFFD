#include <medusa/bits/domains/FindBalancedSupport.hpp>

/**
 * @file
 * Implementation of balanced supports.
 */

namespace mm {

FindBalancedSupport::FindBalancedSupport(int min_support, int max_support) :
        min_support_(min_support), max_support_(max_support), for_which_(), search_among_(),
        force_self_(false) {}

FindBalancedSupport& FindBalancedSupport::forNodes(indexes_t for_which) {
    this->for_which_ = std::move(for_which);
    return *this;
}

FindBalancedSupport& FindBalancedSupport::searchAmong(indexes_t search_among) {
    this->search_among_ = std::move(search_among);
    return *this;
}

FindBalancedSupport& FindBalancedSupport::forceSelf(bool b) {
    force_self_ = b;
    return *this;
}

FindBalancedSupport& FindBalancedSupport::minSupportSize(int size) {
    min_support_ = size; return *this;
}

FindBalancedSupport& FindBalancedSupport::maxSupportSize(int size) {
    max_support_ = size; return *this;
}

}  // namespace mm
