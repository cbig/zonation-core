# OGC ghpages-SWG-template
This is the gh-pages jekyll template for all of the SWG that have an associated website.
To learn how to use this template you can start here: https://jekyllrb.com/docs/home/

## Example build site:
https://opengeospatial.github.io/ghpages-template/

## Making the template your own:
First thing to do is fork this repo to your repo and set it as the gh-pages branch of your repo

You will want to edit the `_config.yml` file with the appropriate settings.
The most obvious change being the title:StandardName. You will want to change this to reflect the name of the standard for which you are working with.

The `_data` directory contains the files that need to edited for each section of the site:

`about.yml`
  This is the section below the header section.
   
`covertopic.yml`
  This is the first section and contains the logos, links and latest news.

`examples.yml`
  This section is below the about section.  

`news.yml`
  This content is found in the _Latest News_ block. The `feed.xml` file lists the items in
  the `news.yml` file.
