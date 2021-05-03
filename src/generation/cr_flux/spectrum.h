#include <TF1.h>

#include <memory>

namespace Utils {
std::unique_ptr<TF1> GetSpectrum(unsigned int Z, unsigned int A);

}