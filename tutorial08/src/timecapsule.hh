template<typename RF>
class TimeCapsule
{
  RF time;
public:
  TimeCapsule ()
  {
    time = 0.0;
  }

  TimeCapsule (RF t)
  {
    time = t;
  }

  // get time
  RF getTime () const
  {
    return time;
  }

  // store time
  void setTime (RF t)
  {
    time = t;
  }
};
