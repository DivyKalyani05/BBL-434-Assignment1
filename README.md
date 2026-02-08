30th Jan
for the main solution i have used a deletion technique where i first found out the ORI via AT rich site method rather than the GC Skew method , as i researched and found out that the GC skew method gives out a weak signal for smaller sequences like plasmids.
I then loaded the markers.tab and deleted the MCS which were absent in design.txt but present in markers.tab file , while protecting the ORI Region... this is similar to the example test case provided
I used this as my primary solution cause when i read about plasmids i found that they have more things in their sequence ,apart from ORI region MCS and antibiotic inhibitor sights, which remain preserved in this solution.

8th Feb 
In my latest commit on 8th Feb i also added the widely used easier method , though wrong in my opinion to my submission as an alternate solution , where I detect the ORI region and extract the replication site , and ultimately add all the required MCS and 
antibiotic resistance sites to make a new plasmid.
This method though will cause the plasmid to not work in that particular microorganism according to my research.
