#ifndef __SIMULATOR_H__
#define __SIMULATOR_H__


#include <random>
#include <thread>
#include <memory>

#include "protocol.h"
#include "record.h"


// TODO move moran over to this system so that genetic algorithms etc
// can be run on either model in the future.
class Simulator {
public:
  virtual std::thread run(Protocol) = 0;
  virtual Record &get_record() = 0;
  // get a shared_ptr to a new copy of the simulator
  virtual std::shared_ptr<Simulator> get_instance() = 0;
private:
};


#endif
