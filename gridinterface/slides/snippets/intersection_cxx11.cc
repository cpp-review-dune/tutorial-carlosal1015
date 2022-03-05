typename SomeGridView::template Codim <0>::Iterator
    it = gridview.template begin<0>();

// iterate over intersections of current entity
for (auto iit = gridview.ibegin(*it);
     iit != gridview.iend(*it); ++iit)
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
