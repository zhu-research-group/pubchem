
from urllib import request, error
import xml.etree.cElementTree as ET


PUBCHEM_BASE = 'http://pubchem.ncbi.nlm.nih.gov/rest/pug/'
NAMESPACE = {'pubchem': 'http://pubchem.ncbi.nlm.nih.gov/pug_rest'}

class AID:
    
    def __init__(self, aid):
        self.aid = str(aid)

    def get_property(self, property):
        """ Returns a particular assay property  """
        try:
            xml = request.urlopen('{0}assay/aid/{1}/summary/XML'.format(PUBCHEM_BASE, self.aid))
        except error.HTTPError as err:
            raise Exception("Sorry, could not pull info for PubChem AID {0}.  Error: {1}".format(self.aid, err))
        except error.URLError as err:
            raise Exception("Sorry, could not pull info for PubChem AID {0}.  Error: {1}".format(self.aid, err))
        except TimeoutError:
            raise Exception("Sorry, could not pull info for PubChem AID {0}.  Error: {1}".format(self.aid, err))
        tree = ET.parse(xml)
        root = tree.getroot()
        attribs = root.findall('pubchem:AssaySummary/pubchem:{0}'.format(property), NAMESPACE)
        property_text = [attrib.text for attrib in attribs if attrib.text]
        return property_text

