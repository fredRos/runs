#include "partitions.h"

#include <iostream>

int main()
{
    for (PartitionGenerator gen(6); gen; ++gen)
        std::cout << *gen << std::endl;

    std::cout << "-----------\n";

    for (KPartitionGenerator gen(6, 3); gen; ++gen)
        std::cout << *gen << std::endl;

    return 0;
}
