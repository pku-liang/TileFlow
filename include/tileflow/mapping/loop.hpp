#pragma once 

#include "mapping/loop.hpp"

namespace loop {

namespace TileFlow{

    class Descriptor: public loop::Descriptor {
    public: 
        std::string name_;
        void Print(std::ostream& out, bool long_form = true) const;
    };

} // namespace TileFlow 

} // namespace loop 

