#ifndef COMMUNICATE_HH
#define COMMUNICATE_HH

template <typename IndexSet, typename Data>
class ExampleDataHandle
  : public Dune::CommDataHandleIF< ExampleDataHandle< IndexSet, Data>, typename Data::value_type >
{
  const IndexSet& iset_;
  Data& data_;
  int cdim_;

public:
  ExampleDataHandle(const IndexSet& iset, Data& data, int cdim) :
    iset_(iset), data_(data), cdim_(cdim)
  {}

  // Returns true if data for this codimension should be communicated
  bool contains (int dim, int codim) const
  {
    // TODO: Add code here
  }

  // Returns true if size per entity of given dim and codim is a constant
  bool fixedSize (int dim, int codim) const
  {
    // TODO: Add code here
  }

  // How many objects of type Data::value_type will we send for this
  // entity?
  template<class EntityType>
  size_t size (EntityType& e) const
  {
    // TODO: Add code here
  }

  // Pack data from user to message buffer. With buff.write(..) you
  // can write something into the buffer and with
  // data_[iset_.index(e)] you can get the entry of the data vector
  // corresponding to e.
  template<class MessageBuffer, class EntityType>
  void gather (MessageBuffer& buff, const EntityType& e) const
  {
    // TODO: Add code here
  }

  // Unpack data from message buffer to user. Hints:
  // - We communicate data of type Data::value_type.
  // - buff.read(x) reads something from the buffer and stores it in x
  // - Acces the data_ vector like above
  template<class MessageBuffer, class EntityType>
  void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
  {
    // TODO: Add code here
  }
};


template <typename GV>
void communicate(const GV& gv, Dune::ParameterTree& ptree){

  // Get a set of indices for entities of current grid view. Note: You
  // only get indices for faces of this processor.
  typedef typename GV::IndexSet IndexSet;
  const IndexSet& indexSet = gv.indexSet();

  // Get the amount of codim 1 entities on this rank and create a data
  // vector of this size. In the 2D case codim 1 entities are the
  // edges.
  int cdim=ptree.get("communication.cdim",(int)1);
  const int dataSize = indexSet.size(cdim);
  std::vector<int> data(dataSize, 0.0);

  // Store the rank of the current processor as data for each egde.
  int myrank = gv.comm().rank();
  for(int i=0 ; i<dataSize; ++i){
    data[i] = myrank;
  }

  // Data handle to communicate the data of each edge to the next processor.
  using DH = ExampleDataHandle<IndexSet, decltype(data)>;
  DH dh(indexSet, data, cdim);

  // Different communication types for DataHandles:
  //
  // InteriorBorder_InteriorBorder_Interface: send/receive interior and border entities
  // InteriorBorder_All_Interface:            send interior and border, receive all entities
  // Overlap_OverlapFront_Interface:          send overlap, receive overlap and front entities
  // Overlap_All_Interface:                   send overlap, receive all entities
  // All_All_Interface:                       send all and receive all entities
  int communicationType = ptree.get("communication.type",(int)5);
  // switch (communicationType){
  // case 1: gv.communicate(dh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication); break;
  // case 2: gv.communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication); break;
  // case 3: gv.communicate(dh,Dune::Overlap_OverlapFront_Interface,Dune::ForwardCommunication); break;
  // case 4: gv.communicate(dh,Dune::Overlap_All_Interface,Dune::ForwardCommunication); break;
  // default: gv.communicate(dh,Dune::All_All_Interface,Dune::ForwardCommunication);
  // }

  // Calculate the sum of the data vector
  int sum(0);
  for(std::size_t i=0 ; i<data.size(); ++i){
    sum += data[i];
  }

  // If we are on rank 0 print the results.
  if (myrank==0){
    std::cout << std::endl;
    std::cout << "== Output for rank " << myrank << std::endl;
    std::cout << std::endl;
    std::cout << "If you run the program with two processes the sum is equal to the number" << std::endl;
    std::cout << "of codim entities where communication happened. The size  of the data" << std::endl;
    std::cout << "vector is equal to the number of all codim entities of this processor." << std::endl;
    std::cout << std::endl;
    std::cout << "Sum of data vector: " << sum << std::endl;
    std::cout << "Size of data vector: " << data.size() << std::endl;
  }
}

#endif
