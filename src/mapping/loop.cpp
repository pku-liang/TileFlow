#include "tileflow/mapping/loop.hpp"

namespace loop{

namespace TileFlow {

void Descriptor::Print(std::ostream& out, bool long_form) const
{
  if (long_form)
  {
    out << "for " << name_ << " in [" << start << ":" << end;
    if (residual_end != end)
      out << "," << residual_end;
    out << ")";
    if (loop::IsSpatial(spacetime_dimension))
    {
      if (loop::IsSpatialX(spacetime_dimension))
        out << " (Spatial-X)";
      else
        out << " (Spatial-Y)";
    }
  }
  else
  {
    out << "(" << name_ << "," << end;
    if (residual_end != end)
      out << "," << residual_end;
    if (loop::IsSpatial(spacetime_dimension))
    {
      if (loop::IsSpatialX(spacetime_dimension))
        out << ",spX";
      else
        out << ",spY";
    }
    out << ") ";
  }
}

} // namespace TileFlow 

} // namespace loop