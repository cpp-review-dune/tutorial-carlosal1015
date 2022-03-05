typename SomeGridView::template Codim <0>::Iterator
    it = gridview.template begin<0>();
typedef typename SomeGridView::IntersectionIterator
    IntersectionIterator;

// iterate over intersection of current entity
IntersectionIterator iitend = gridview.iend(*it);
for (IntersectionIterator iit = gridview.ibegin(*it);
     iit != iitend; ++iit)
{
   // neighbor intersection
   if (iit->neighbor())
   {
       // do something ...
   }
   // boundary intersection
   if (iit->boundary())
   {
       // do something else ...
   }
}
