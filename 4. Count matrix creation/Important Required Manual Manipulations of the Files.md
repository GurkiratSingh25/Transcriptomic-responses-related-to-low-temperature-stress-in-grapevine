#Originally written by: Cullen W. Dixon

Step 1 - 
Manually add in " gene_biotype vitis_labrusca;" to the end of column 9 lines which contain the word "gene" in the 3rd column.
You can put in any genus_species.
I recommend doing this in excel using filters and autofill.
You can probably figure out a way to do this using the help of ChatGPT now but this was before that was a tool and finding something in one column and then acting
on another is not an inconsequential coding problem to resolve.

Step 2 - 
Manually remove all instances of doubled-up quotation marks.
Example = ""text is here"" 
Without removal of these the program gets really confused.
I recommend doing this in sublime using find and replace in which you'll search for {""} and replace with {"}.

Step 3 - 
Change line endings in sublime to 'Unix' from 'Windows' since we have (presumably) made changes using excel (unless you did it some other way).

Step 4 -
Bask in the glory of having circumvented a major roadblock in the pipeline.
