#pragma once 

#include "mapping/loop.hpp"

namespace loop {

namespace TimeloopX{

    class Descriptor: public loop::Descriptor {
    public: 
        std::string name_;
        void Print(std::ostream& out, bool long_form = true) const;
    };

} // namespace TimeloopX 

} // namespace loop 

