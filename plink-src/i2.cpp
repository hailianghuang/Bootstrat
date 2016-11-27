

//////////////////////////////////////////////////////////////////
//                                                              //
//           PLINK (c) 2005-2009 Shaun Purcell                  //
//                                                              //
// This file is distributed under the GNU General Public        //
// License, Version 2.  Please see the file COPYING for more    //
// details                                                      //
//                                                              //
//////////////////////////////////////////////////////////////////


// TODO:
//  1. Implement equivalence classes
//  2. Implement and check behavior of joint fields
//     -- how does this interact with lookup and replace options?
//  3. Add match command


#include "idhelp.h"
#include "options.h"
#include "helper.h"
#include "nlist.h"
#include <iomanip>

class IDFile;

class IDField
{
public:

  string name;

  bool null;
  bool attribute;
  bool equiv;
  bool joint;


  // Equivalence IDs
  map<string,string > eqid;
  set<string> aliasList;
  set<string> masterList;

  
  IDField()
  {
    attribute = false;
    null = false; 
    equiv = false;
    joint = false;
  }
  
  bool operator< (const IDField & b) const
    {
      if ( name < b.name ) 
	return true;
      return false;
    }

  bool operator== (const IDField & b) const
    {
      return name == b.name; 
    }
};

class IDFile 
{
public:

  string filename;

  int uniqFieldCount;
  
  bool hasHeader;
  
  set<string> missingValues;
  
  vector<IDField*> fields;
  vector< vector<int> > joint;
  vector< vector<int> > equiv;
  map<string,string> injections;

  IDFile() 
    { 
      hasHeader = false;
      missingValues.insert(".");
      uniqFieldCount = 0;
    }
};

class IDValue
{
public:
  
  IDValue()
  {
    value = "";
    field = NULL;
    jointValue = "";
  }

  void updateAlias()
  {
    map<string,string>::iterator f = field->eqid.find( value );
    if ( f != field->eqid.end() )
      value = f->second;
  }

  IDField * field;
  
  string jointValue;
  
  // This is main default value
  string value;

  
  bool operator< (const IDValue & b) const
  {
    if ( field->name < b.field->name )
      return true;
    else if ( field->name > b.field->name )
      return false;
    if ( field->joint )
      return jointValue < b.jointValue;
    else
      return value < b.value; 

  }  


  bool operator!= (const IDValue & b) const
  {
    return ! ( *this == b );
  }


  bool operator== (const IDValue & b) const
  {
        
    // These must refer to the same field to be 
    // comparable

    if ( field->name != b.field->name ) 
      return false;

    if ( field->joint )
      {
	return jointValue == b.jointValue;
      }
    else
      {
	return value == b.value;
      }
  }


  // NOT NEEDED  
//   bool inconsistent(IDValue * b ) 
//   {
//     if ( field->name == b->field->name && 
// 	 value != b->value )
//       return true;
//     return false;    
//   }

//   bool consistent(IDValue * b ) 
//   {
//     if ( field->name == b->field->name && 
// 	 value == b->value )
//       return true;
//     return false;    
//   }
  
//   bool sameField(IDValue * b ) 
//   {
//     return field->name == b->field->name;
//   }

};


class IDGroup
{
public:

  vector<IDValue*> values;
  IDFile * file;

  bool resolved;

  IDGroup()
  {
    resolved = false;
  }

  void display()
  {
    cout << "File = " << file->filename 
	 << ", resolved = " 
	 << resolved << "\n";
    for (int k=0; k<values.size(); k++)
      cout << "\t" << values[k]->field->name 
	   << " = " << values[k]->value 
	   << "\n";
    cout << "\n";
  }

};




// class UniquePerson
// {
// public:
//   set<IDValue> id;
//   bool matches(IDGroup &);
// };


// // Map one or more people given a set of IDs
// set<UniquePerson*> findPerson( IDGroup & id , vector<UniquePerson*> & people )
// {
//   set<UniquePerson*> matches;
//   for ( int i = 0 ; i < people.size() ; i++ )
//     {
//       if ( people[i]->matches( id ) )
// 	matches.insert( people[i] );
//     }
//   return matches;
// } 

// bool UniquePerson::matches( IDGroup & id )
// {
//   // Does this person match with the ID value?
//   return false;
// }


void Plink::idHelp()
{
  
  // All people -- collate harmonized IDs in here
  //  vector<UniquePerson> people;
  
  // Files containing IDs
  vector<IDFile> files;
  
  // The actual IDs we are matching on
  set<IDField> fields;
  map<string,IDField*> fieldMap;
  set<IDField>::iterator iField;
  
  // The basic data we read in, then try to resolve into UniquePerson's
  vector<IDGroup> idgroup;
    
  // The lookup table
  map<IDValue,set<IDGroup*> > idmap;
   

  // 1. Read in dictionary, and make the fields and files

  // Contains files (and can be full path) and description of each field
  // {file name}   {col names } : { rules } 
  // e.g.

  // ../files/id1.txt FID IID FID1 FID2 : uniq=FID,IID uniq=FID1,IID2 
  // ../names/id.lst ID2 ID23 : missing=NA,---,0 
  // ../names/id2.lst ID3 : equiv 
  // 

  // note: for "equiv" files, assume all IDs on same line are equivalent, 
  // only one ID can be specified here

  
  printLOG("ID helper, with dictionary [ " + par::idhelp_dictionary + " ]\n"); 
  checkFileExists( par::idhelp_dictionary );  

  map<string,IDFile> dict;
 
  set< set<string> > uniqueFields;
  map<string,set<string> > jointMap;

  set<IDField*> isJoint;
  vector< set<IDField*> > jointField;

  set<string> attribFields;
  
  ifstream DICT( par::idhelp_dictionary.c_str() , ios::in );

  while ( ! DICT.eof() )
    {
      vector<string> tokens = tokenizeLine( DICT );
      // Needs atleast three fields:
      // filename, and two IDs
      
      if ( tokens.size() == 0 ) 
	continue;
      
      if ( tokens.size() < 3 || tokens[1] == ":" || tokens[2] == ":" ) 
	error("Expecting at least 3 fields w/ 2 IDs in every dictionary row\n");
      
      IDFile d;
      checkFileExists( tokens[0] );
      d.filename = tokens[0];
      
      int p = 1;

      while (1)
	{
	  
	  if ( tokens[p] == ":" )
	    {
	      ++p;
	      break;
	    }
	  
	  // Is this a new field name?
	  IDField f;
	  f.name = tokens[p];
	  if ( f.name == "." ) 
	    f.null = true;
	  iField = fields.find( f );
	  
	  if ( iField == fields.end() )
	    {
	      fields.insert( f );
	      iField = fields.find( f );
	    }
	  
	  // Track that this field was found in this file
	  
	  IDField * ip = (IDField*)&( *iField );
	  d.fields.push_back( ip );
	  ++d.uniqFieldCount;
 
	  // Consider the next field in this file
	  if ( ++p == tokens.size() )
	    break;
	}
      

      // Now read any rules 
      if ( p < tokens.size() )
	{
	  while (1)
	    {
	      string cmd = tokens[p];

	      // Commands
	      // attrib=X,Y,Z
	      // joint=X,Y
	      // header
	      // missing=NA,-9

	      if ( cmd.size() > 6 && cmd.substr(0,6) == "joint=" )
		{
		  string u = cmd.substr(6);
		  // This should contain at least two ID fields
		  // These fields must always then appear together in 
		  // any file that features at least one
		  
		  NList tlist(0);
		  vector<string> ids = tlist.deparseStringList( u );
		  if ( ids.size() < 2 )
		    error("Problem with specification of : " + cmd );
		  set<string> t;
		  for (int i = 0 ; i < ids.size() ; i++)
		    {
		      t.insert(ids[i]);
		    }
		  
		  for (int i = 0 ; i < ids.size() ; i++)
		    jointMap.insert( make_pair( ids[i] , t ) ); 
				      
		  set<IDField*> jointf;

		  for (int i = 0 ; i < ids.size() ; i++)
		    {
		      IDField f;
		      f.name = ids[i];
		      iField = fields.find( f );
		      if ( iField == fields.end() )
			{
			  error("Could not find field " + ids[i] + " which is specified as joint");
			}
		      
		      IDField * ip = (IDField*)&( *iField );
		      ip->joint = true;
		      jointf.insert( ip );
		      isJoint.insert( ip );
		    }
		  jointField.push_back( jointf );

		}

	      
	      if ( cmd.size() > 7 && cmd.substr(0,7) == "attrib=" )
		{
		  string u = cmd.substr(7);
		  NList tlist(0);
		  vector<string> ids = tlist.deparseStringList( u );
		  for (int i = 0 ; i < ids.size() ; i++)  
		    attribFields.insert( ids[i] );
		}
	      
	      if ( cmd.size() > 8 && cmd.substr(0,8) == "missing=" )
		{
		  string u = cmd.substr(8);
		  NList tlist(0);
		  vector<string> ids = tlist.deparseStringList( u );
		  for (int i = 0 ; i < ids.size() ; i++)  
		    d.missingValues.insert( ids[i] );
		}
	      
	      if ( cmd.size() > 4 && cmd.substr(0,4) == "set:" )
		{
		  string u = cmd.substr(4);
		  bool okay = true;

		  if ( u.find("=") == string::npos )
		    okay = false;
		  else
		    {
		      string u1 = u.substr(0,u.find("="));
		      string u2 = u.substr(u.find("=")+1);
		      if ( u1.size() < 1 )
			okay = false;
		      if ( u2.size() < 1 ) 
			okay = false;
		      
		      if ( okay ) 
			{
			  d.injections.insert(make_pair(u1,u2));
			  IDField f;
			  f.name = u1;
			  if ( f.name == "." )
			    error("Cannot set field value name to .");
			  
			  iField = fields.find( f );
			  
			  if ( iField == fields.end() )
			    {
			      fields.insert( f );
			      iField = fields.find( f );
			    }
			  
			  // Track that this field was found in this file
			  IDField * ip = (IDField*)&( *iField );
			  d.fields.push_back( ip );
			  ++d.uniqFieldCount;
			  
			}
		      
		    }

		  if ( ! okay ) 
		    error("Badly formed set:X=V command");
		  
		}

	      if( cmd == "header" || cmd == "hasHeader" )
		d.hasHeader = true;
	      
	      if ( ++p == tokens.size() )
		break;
	    }
	}
      
      files.push_back(d);
      // Next line 
    }
  
  DICT.close();  

  printLOG("Read " + int2str( fields.size() ) + " unique fields\n");
  
  set<IDField>::iterator i = fields.begin();
  while ( i != fields.end() )
    {
      fieldMap.insert(make_pair( i->name , (IDField*)&(*i) ));
      ++i;
    }


  /////////////////////////////////////////////////
  // If set fields, check either joint or attrib
  
  for ( int f = 0 ; f < files.size() ; f++ ) 
    {
      map<string,string>::iterator i = files[f].injections.begin();
      while ( i != files[f].injections.end() )
	{
	  IDField * f = fieldMap.find( i->first )->second;
	  if ( ! ( f->joint || f->attribute ) )
	    error("Any set:field=value should be an attribute or a joint field");
	  ++i;
	}
    }

  ////////////////////////////////////////////////
  // If joint fields, check always all specified

  for (int j = 0 ; j < jointField.size(); j++ )
    {
      set<IDField*> & jf = jointField[j];
      for ( int f = 0 ; f < files.size() ; f++ ) 
	{
	  for ( int j = 0; j < files[f].fields.size(); j++)
	    {
	      if ( jointMap.find( files[f].fields[j]->name ) != jointMap.end() )
		{
		  map<string,set<string> >::iterator i = jointMap.find( files[f].fields[j]->name );
		  set<string> & ss = i->second;
		  set<string>::iterator is = ss.begin();
		  while ( is != ss.end() )
		    {
		      bool okay = false;
		      for ( int j = 0; j < files[f].fields.size(); j++)
			if ( *is == files[f].fields[j]->name ) 
			  okay =true;
		      if ( ! okay ) 
			error("Need to specify all joint fields here");
		      ++is;
		    }
		}
	    }
	}
    }


  /////////////////////////////////
  // 2. Read in and index all IDs

  for ( int f = 0 ; f < files.size() ; f++ ) 
    {


      IDFile * file = &files[f];

      printLOG("Reading [ " + files[f].filename + " ] with fields : ");

      for ( int j = 0 ; j < files[f].fields.size(); j++ )
	{
	  if (j>0 )
	    printLOG(", ");
	  printLOG( files[f].fields[j]->name );
	}


      ////////////////////////////
      // Find any equivalence sets

      bool foundEquiv = false;
      map<string, vector<int> > equiv;
      map<string, map<int,int> > equivMap;
      for ( int j = 0 ; j < files[f].fields.size(); j++)
	{
	  string n = files[f].fields[j]->name;
	  map<string, vector<int> >::iterator i = equiv.find( n );
	  if ( i == equiv.end() )
	    {
	      vector<int> t;
	      t.push_back(j);
	      equiv.insert(make_pair( n , t ));
	    }
	  else
	    {
	      i->second.push_back(j);
	      files[f].fields[j]->equiv = true;
	      foundEquiv = true;	      
	    }
	}
      if ( foundEquiv ) 
	{
	  printLOG(" : ");
	  map<string, vector<int> >::iterator i = equiv.begin();
	  while ( i != equiv.end() )
	    {
	      if ( i->second.size() > 1 )
		{
		  map<int,int> tmap;
		  
		  printLOG(" "+files[f].fields[i->second[0]]->name + "("); 
		  for (int k = 0 ; k < i->second.size(); k++)
		    {
		      if ( k>0 ) 
			{
			  if ( k==1 ) 
			    printLOG("<-");
			  else
			    printLOG(",");
			  tmap.insert(make_pair( i->second[k] , i->second[0] ) );
			}
		      printLOG(int2str( i->second[k]+1 ));
		    }
		  printLOG(")");
		  
		  equivMap.insert( make_pair( files[f].fields[i->second[0]]->name , tmap ));
		}
	      
	      ++i;
	    }
	}
    
      
      printLOG("\n");

      
            

      ////////////////////////////////////
      // Read the raw data
      

      ifstream ID1( files[f].filename.c_str() , ios::in );
      
      if ( files[f].hasHeader )
	{
	  vector<string> header = tokenizeLine( ID1 );
	}

      while ( !ID1.eof() )
	{

	  vector<string> tokens = tokenizeLine( ID1 );
	  if ( tokens.size() == 0 ) 
	    continue;

	  // Insert SET:X=Y values here
	  map<string,string>::iterator i = files[f].injections.begin();
	  while ( i != files[f].injections.end() )
	    {
	      tokens.push_back( i->second );
	      // Note -- this won't be same oreder -- need to check/fix this
	      ++i;

	    }

	  if ( tokens.size() != file->uniqFieldCount )
	    {
	      printLOG("\n\nIn [ " + file->filename 
		       + " ] encountered a row with the wrong number of fields\n");
	      printLOG("Found " + int2str( tokens.size() ) 
		       + " fields but expecting " + int2str( file->fields.size() ) + "\n");
	      int mx = tokens.size() > file->uniqFieldCount ? tokens.size() : file->uniqFieldCount ;

	      for (int j = 0 ; j < mx ; j++)
		{
		  if ( j < file->uniqFieldCount )
		    printLOG( "   " + file->fields[j]->name + ":\t" );
		  else
		    printLOG( "   {?}:\t" );
		  
		  if ( j < tokens.size() )
		    printLOG( "   " + tokens[j] + "\n" );
		  else
		    printLOG( "   {?}\n" );

		}
	       

	      error("Problem with [ " + file->filename + " ]\n" );
	    }

	  IDGroup g;

	  // Track which file this group of IDs came form
	  g.file = file;

	  // Track what the original eq-value is for this line
	  map<IDField*,string> originalEquivalence;
	  
	  for ( int j = 0 ; j < files[f].fields.size(); j++ )
	    {
	      
	      // Is this a missing value?
	      
	      if ( files[f].missingValues.find( tokens[j] ) != files[f].missingValues.end() )
		continue;
	      	      
	      IDField * myField = files[f].fields[j];
	      
	      // Is this field an equivalence set? 

	      bool needToStore = true;

	      if ( myField->equiv )
		{
		  map<IDField*,string>::iterator i = originalEquivalence.find( myField );
		  if ( i == originalEquivalence.end() )
		    {
		      
		      set<string>::iterator k = myField->aliasList.find( tokens[j] );
 		      if ( k != myField->aliasList.end() )
 			error("Alias specified more than once in " 
			      + files[f].filename + " : " 
			      + myField->name + " " + tokens[j]);
		      myField->masterList.insert( tokens[j] );
		      originalEquivalence.insert(make_pair( myField, tokens[j] ));		      
		    }
		  else
		    {
		      
		      // We need to add an equivalence value to this field?
		      map<string, string>::iterator k = myField->eqid.find( tokens[j] );
		      
		      // We only let each alias be specified once
		      if ( k != myField->eqid.end() )
			error("Alias specified more than once in " 
			      + files[f].filename + " : " 
			      + myField->name + " " + tokens[j]);
		      
		      k = myField->eqid.find( i->second );
 		      if ( k != myField->eqid.end() )
 			error("Alias specified more than once in " 
			      + files[f].filename + " : " 
			      + myField->name + " " + i->second);
		      
		      if ( myField->masterList.find( tokens[j] ) != myField->masterList.end() )
			error("Alias specified more than once in " 
			      + files[f].filename + " : " 
			      + myField->name + " " + tokens[j]);
		      
		      
		      // Keep track of this alias
		      myField->eqid.insert( make_pair ( tokens[j] , i->second ) );
		      
		      myField->aliasList.insert( tokens[j]  );
		      
		      // And now we do not need to store this particular value
		      // as a distinct field
		      
		      needToStore = false;
		    }

		}
	      
	      if ( needToStore )
		{
		  IDValue * v = new IDValue;
		  v->field = myField;		  
		  v->value = tokens[j];	      
		  g.values.push_back(v);
		}
	      
	    }

	  idgroup.push_back(g);
	}
      
      ID1.close();
    }

  
  /////////////////////////////////
  // Done reading in the raw data

  if ( attribFields.size() > 0 )
    {
      printLOG("   Attribute fields: ");
      set<string>::iterator i = attribFields.begin();
      while ( i != attribFields.end() )
	{
	  if ( fieldMap.find( *i ) == fieldMap.end() )
	    error("Cannot find specified attribute " 
		  + *i 
		  + " -- please check your dictionary file");

	  // Set as an attribute field
	  fieldMap.find(*i)->second->attribute = true;

	  printLOG( *i + " " );
	  ++i;
	}
      printLOG("\n");
    }
  

  if ( uniqueFields.size() > 0 )
    {
      string problem = "";
      printLOG("   Joint fields: ");
      set<set<string> >::iterator i = uniqueFields.begin();
      while ( i != uniqueFields.end() )
	{
	  printLOG("{ ");
	  set<string>::iterator j = i->begin();
	  while ( j != i->end() )
	    {

	      if ( attribFields.find( *j ) != attribFields.end() )
		problem += *j + " cannot be both an attribute and a joint ID\n";

	      printLOG( *j + " " );
	      ++j;
	    }
	  printLOG("} ");
	  ++i;
	}
      printLOG("\n");

      if ( problem.size() > 0  ) 
	error( problem );
    }



 
  //////////////////////////////////////////////////////////////
  // 1. Swap in preferred values for any equiv fields
  
  for ( int g = 0 ; g < idgroup.size(); g++ )
    {
      IDGroup * group = &idgroup[g];
      for (int j = 0 ; j < group->values.size(); j++)
	if ( group->values[j]->field->equiv )
	  group->values[j]->updateAlias();
    }



  //////////////////////////////////////////////////////////////
  // 1b. Compile joint fields into joint values

  for (int j = 0 ; j < jointField.size(); j++ )
    {
      set<IDField*> & jf = jointField[j];
      
      // Does this joint field exist witin the 
      // group?
      
      for ( int g = 0 ; g < idgroup.size(); g++ )
	{
	  IDGroup * group = &idgroup[g];
	  
	  // Does this group contain a joint field?

	  bool hasJoint = false;
	  map<IDField*,int> mapback;
	  for (int j = 0 ; j < group->values.size(); j++)
	    {
	      if ( jf.find( group->values[j]->field ) != jf.end() )
		{
		  hasJoint = true;
		  mapback.insert( make_pair( group->values[j]->field , j ) );
		}
	    }
	  
	  if ( ! hasJoint ) 
	    continue;


	  string jointValue = "";
	  bool doneFirst = false;
	  bool jointMissing = false;
	  
	  set<IDField*>::iterator k = jf.begin();
	  while ( k != jf.end() )
	    {
	      
	      map<IDField*,int>::iterator mi = mapback.find( *k );
	      if ( mi == mapback.end() )
		jointMissing = true;
	      else
		{
		  if ( doneFirst ) 
		    jointValue += "+" + group->values[ mi->second ]->value;
		  else
		    {
		      jointValue += group->values[mi->second]->value;
		      doneFirst = true;
		    }
		}
	      ++k;
	    }
	  

	  
	  // Update values accordingly
	  vector<bool> mask( group->values.size(), false);

	  set<IDField*>::iterator k2 = jf.begin();
	  while ( k2 != jf.end() )
	    {
	      
	      map<IDField*,int>::iterator mi = mapback.find( *k2 );
	      if ( mi != mapback.end() )
		{
		  int j = mi->second;
		  
		  if ( jointMissing ) 
		    mask[j] = true;
		  else
		    group->values[j]->jointValue = jointValue;
		}
	      ++k2;
	    }
	  if ( jointMissing )
	    {
	      vector<IDValue*> newValues = group->values;
	      group->values.clear();
	      for ( int i = 0 ; i < mask.size() ; i++)
		{
		  if ( ! mask[i] ) 
		    group->values.push_back( newValues[i] );
		}

	    }
	}  

    }



  ////////////////////////////////
  // 2. Create the idmap


  for ( int g = 0 ; g < idgroup.size(); g++ )
    {
      IDGroup * group = &idgroup[g];
      for (int j = 0 ; j < group->values.size(); j++)
	{

	  IDValue & v = *(group->values[j]);
	  map<IDValue,set<IDGroup*> >::iterator i = idmap.find( *(group->values[j]) );
	  if ( i == idmap.end() )
	    {
	      set<IDGroup*> t;
	      t.insert(group);
	      idmap.insert(make_pair( *(group->values[j]) , t) );
	    }
	  else
	    {
	      i->second.insert( group );
	    }
	}
    }



  ///////////////////////////////////////////////
  // 3. Attempt to resolve into a single table
  
  
  bool okay = true;
  map<string,string> problem ;

  while (1)
    {

      bool allDone = true;

      for ( int g = 0 ; g < idgroup.size(); g++ )
	{

	  IDGroup * group = &idgroup[g];
	  
	  // Has this group already been assigned to a person?
	  
	  if ( group->resolved )
	    continue;
	  
	  
	  // Find all other groups (resolved or otherwise) that this group 
	  // matches with, but ignoring attributes
	  
	  set<IDGroup*> matches;
	  
	  for (int j = 0 ; j < group->values.size(); j++)
	    {
	      // Skip matching on attributes
	      if ( group->values[j]->field->attribute )
		continue;

	      map<IDValue, set<IDGroup*> >::iterator i = idmap.find( *(group->values[j]) );
	      if ( i != idmap.end() )
		{
		  set<IDGroup*>::iterator i2 = i->second.begin();
		  while ( i2 != i->second.end() )
		    {
		      matches.insert( *i2 );
		      ++i2;
		    }
		}
	    }
	  

	  //////////////////////////////
	  // Merge into the key group
      	  
	  // Make a set of the key groups IDValues

	  map<IDField*,IDValue*> keyValues;
	  for ( int j = 0; j < group->values.size(); j++ )
	    keyValues.insert( make_pair( group->values[j]->field, group->values[j] ));
	  
	  set<IDGroup*>::iterator i0 = matches.begin();
	  while ( i0 != matches.end() )
	    {
	      
	      if ( *i0 == group ) 
		{
		  ++i0;
		  continue;
		}
	      
	      
	      // Step through all the values in this matching group
	      
	      for ( int k = 0; k < (*i0)->values.size(); k++)
		{
		  // Does the key have this field? 

		  IDField * f = (*i0)->values[k]->field;

		  map<IDField*,IDValue*>::iterator i = keyValues.find( f );
		  
		  if ( i == keyValues.end() )
		    {
		      // Insert this key value into place
		      IDValue * t = new IDValue;
		      t->field = f;
		      t->value = (*i0)->values[k]->value;
		      group->values.push_back(t);
		      allDone = false;
		    }
		  else
		    {

		      // Something already exists -- check it is not inconsistent
		      // as if is an ID, i.e. attributes are allowed to not be 
		      // unique by definition
		      
		      if (  *((*i0)->values[k]) != *(i->second) )
			{
			  
			  okay = false;

			  string title = "Two entries are unique [ " + (*i0)->values[k]->field->name + " = ";
			  
			  if ( (*i0)->values[k]->value <=  i->second->value )
			    title += (*i0)->values[k]->value + " and " + i->second->value;
			  else
			    title += i->second->value + " and " + (*i0)->values[k]->value;
			  
			  title += " ] but match on others fields";
			  
			  if ( problem.find( title ) == problem.end() )
			    {
			      string p = "\n a) ";
			      for (int z=0; z<group->values.size(); z++) 
				p += group->values[z]->field->name + "=" 
				  + group->values[z]->value + " ";
			      p += "\n";
			      
			      p += " b) ";
			      for (int z=0; z<(*i0)->values.size(); z++) 
				p += (*i0)->values[z]->field->name + "=" 
				  + (*i0)->values[z]->value + " ";
			      p += "\n\n";
			      
			      problem.insert(make_pair(title,p));
			    }
			  
			}
		    }
		}


	      // We are now done with this IDGroup

	      (*i0)->resolved = true;
	      //(*i0)->display();
	      
	      ++i0;
	    }
	}
      
      if ( allDone )
	break;

    } 

  if ( ! okay ) 
    {
      printLOG("\n\n*** Problems were detected in the ID lists:\n\n"); 
      map<string,string>::iterator p = problem.begin();
      while ( p != problem.end() )
	{
	  printLOG( p->first + p->second );
	  ++p;	
	}
      error("You need to fix the above problems");
    }



  
  
  //////////////////////////////////////////////////////////////
  // 4. Figure out actual transformation required, and perform

  // Rules:  ID cannot contain whitespace or commas
  // Can contain "-", "_", "+", "=", etc.
  // IDs are case-sensitive
  // Values if "." are taken to mean not known, n/a

  // Functions: dump all fields on this person, or group
  if ( par::idhelp_list_aliases )

  // Lookup a single person       --id-lookup ID=27364883-1
  //  or match on an attribute    --id-lookup SITE=Boston
  
  // Dump the entire table        (DEFAULT)
  //     or subset of cols        --id-table ID,CLIN_ID,BSP_ID
  
  // Take an existing file and replace a field  

  //                              --id-replace [header|noheader,field=N,skip|miss|warn|list] mydata.txt ID1 ID2

  // Default behavior is to put 'missing' as the field
  // default = to autodetect a header field
  //           skip  (do not print these lines)
  //           miss  (print a missing ID code)
  //           warn  (do not allow if >1 missing)
  //           list  (print only these lines, ID in file not in DB)

  
  // Take an existing file and replace a field  {file} {field} {old ID} {new ID} 
  //                              --id-replace mydata.txt ID1 CLIN_ID

  // par::idhelp_command = 
  // dump_table,  dump_subtable, replace, lookup

  if ( par::idhelp_list_aliases )
    {
      printLOG("Listing ID equivalents/aliases to [ " + par::output_file_name + ".id.eq ]\n");
      ofstream O1( (par::output_file_name + ".id.eq").c_str() , ios::out );
      
      O1 << setw(20) << "FIELD" << " " 
	 << setw(20) << "PREF" << " "
	 << setw(20) << "EQUIV" << "\n";
      
      map<string,IDField*>::iterator i1 = fieldMap.begin();
      while ( i1 != fieldMap.end() )
	{
	  if ( i1->second->equiv ) 
	    {
	      IDField * f = i1->second;
	      map<string,string>::iterator j = f->eqid.begin();
	      while ( j != f->eqid.end() )
		{
		  O1 << setw(20) << i1->first << " "
		     << setw(20) << j->second << " " 
		     << setw(20) << j->first << "\n";
		  ++j;
		}
	    }
	  ++i1;
	}
      O1.close();
      
      return;
    } 
  
  
  if ( par::idhelp_replace )
    {
      

      
      NList nl(0);

      vector<string> tok = nl.deparseStringList( par::idhelp_replace_string );
      
      if ( tok.size() != 3 ) 
	error("Hmm... internal problem\n");

      checkFileExists( tok[0] );

      if ( fieldMap.find( tok[1] ) == fieldMap.end() ) 
	error("Cannot find " + tok[1] );
      if ( fieldMap.find( tok[2] ) == fieldMap.end() ) 
	error("Cannot find " + tok[2] );
      printLOG("Replacing " + tok[1] + " with " 
	       + tok[2] + " from [ " + tok[0] + " ]\n");
      printLOG("Writing new file to [ " + par::output_file_name + ".rep ]\n");

      bool skipMode = par::idhelp_replace_options.isSet("skip");
      bool missMode = par::idhelp_replace_options.isSet("miss");
      bool warnMode = par::idhelp_replace_options.isSet("warn");
      bool listMode = par::idhelp_replace_options.isSet("list");
      int c = 0;

      if ( skipMode ) 
	{
	  printLOG("Set to skip unmatched observations\n");
	  ++c;
	}

      if ( missMode ) 
	{
	  printLOG("Set to set unmatched observations to NA\n");
	  ++c;
	}

      if ( warnMode ) 
	{
	  printLOG("Set to give error for first unmatched observation\n");
	  ++c;
	}

      if ( listMode ) 
	{
	  printLOG("Set to list only IDs in file but not in database\n");
	  ++c;
	}

      if ( c == 0 )
	printLOG("Set to  keep original value for unmatched observations\n"); 

      if ( c>1 ) 
	error("Can only specify one of [miss|warn|skip|list] options in --id-replace");


      // Do we have a header row?
      bool header = false;
      if ( par::idhelp_replace_options.isSet("header") )
	header = true;

      int rep_field = -1;
      
      
      // If not, we need a number specified
      if ( ! header )
	{
	  string field_str;
	  if ( par::idhelp_replace_options.isSet("field") )
	    field_str = par::idhelp_replace_options.getValue("field");
	  else
	    error("Need to specify field={N} if no header");
	  if ( ! from_string<int>( rep_field, field_str , std::dec ) )
	    error("Problem with field specified in --id-replace options");
	  // Make 0-based
	  --rep_field;
	}
      else
	{
	  // Lookup in header row
	  ifstream IN1( tok[0].c_str() , ios::in );
	  vector<string> tokens = tokenizeLine( IN1 );
	  for (int i = 0 ; i < tokens.size(); i++)
	    {
	      if ( tokens[i] == tok[1] )
		rep_field = i;
	    }
	  IN1.close();
	  IN1.clear();
	  if ( rep_field == -1 ) 
	    error("Could not find " + tok[1] + " in [ " + tok[0] + " ]\n");
	}
      
      // Load input file

      ifstream IN1( tok[0].c_str() , ios::in );
      ofstream OUT1( (par::output_file_name+".rep").c_str() , ios::out );
      int notFound = 0;
      bool readHeader = false;
      if ( ! header ) 
	readHeader = true;
      
      while( ! IN1.eof() )
	{
	  vector<string> tokens = tokenizeLine( IN1 );
	  if ( tokens.size() == 0 )
	    continue;
	  if ( tokens.size() <= rep_field )
	    error("Not enough columns here");
	  
	  bool changed = false;
	  
	  // Deal with header row?
	  if ( ! readHeader ) 
	    {
	      tokens[ rep_field ] = tok[2];
	      readHeader = true;
	      changed = true;
	    }
	  else
	    {
	      
	      // Swap value
	      // From idmap, there should be a single unresovled IDGroup
	      
	      IDValue findField;
	      findField.field = fieldMap.find( tok[1] )->second;
	      findField.value = tokens[ rep_field ];
	      if ( findField.field->equiv )
		findField.updateAlias();
	      	      
	      
	      map<IDValue,set<IDGroup*> >::iterator i = idmap.find( findField );
	      
	      if ( i != idmap.end() )
		{	      
		  int unres = 0;
		  IDGroup * thisGroup = 0;
		  set<IDGroup*>::iterator j = i->second.begin();
		  while ( j != i->second.end() )
		    {
		      if ( ! (*j)->resolved ) 
			{
			  ++unres;
			  thisGroup = *j;
			}
		      ++j;
		    }
		  
		  if ( unres != 1 ) 
		    cout << "Hmm... ******** a problem, unres = " << unres << "\n";
		  
		  for (int i = 0 ; i < thisGroup->values.size(); i++)
		    if ( thisGroup->values[i]->field->name == tok[2] )
		      {
			tokens[ rep_field ] = thisGroup->values[i]->value;
			changed = true;
		      }
		}
	      else
		++notFound;
	    }

	  if ( ! changed ) 
	    {
	      if ( listMode ) 
		{
		  OUT1 << tokens[ rep_field ] << "\n";
		  continue;
		}
	      
	      if ( skipMode )
		continue;
	      
	      if ( warnMode )
		error("Could not find replacement for " + tokens[rep_field] + "\n");
	      
	      if ( missMode )
		tokens[rep_field] = "NA";
	    }
	  
	  if ( listMode ) 
	    continue;

	  // Print out line, which may or may not be modified
	  for (int i = 0 ; i < tokens.size() ; i++ ) 
	    OUT1 << tokens[i] << " ";	  
	  OUT1 << "\n";
	      
	}
     
      if ( notFound > 0 ) 
	printLOG("Could not find matches for " + int2str( notFound ) + " lines\n");

      OUT1.close();
      IN1.close();
      
      return;
    }
  
  
  ///////////////////////////////////////////////////////////////
  // Line up 1 or more files based on the first file
  
  if ( par::idhelp_match )
    {
      
      // e.g. --id-match myfile.fam FID,IID 1,2 file1.txt CLIN_ID 1 file2.txt  
      // in form {file} {id} {col} 
      // If joint IDs, then all must be specified.
      // Cannot specify more than 1 non-joint ID though
      // Can be different IDs in different files
     
      error("Not yet implemented");
      
 
    }
  

  map<IDField*, set<IDValue> > lookupValues;
  set<string>  subsetFields;
  
  if ( par::idhelp_subset )
    {
      NList tlist(0);
      vector<string> ids = tlist.deparseStringList( par::idhelp_subset_string );
      for (int i=0; i<ids.size(); i++)
	{
	  if ( fieldMap.find(ids[i]) == fieldMap.end() )
	    error("Cannot find field " + ids[i] );
	  subsetFields.insert(ids[i]);
	}
    }


  if ( par::idhelp_lookup )
    {

      // Input should be in form of a comma delimited list ID=value, 
      // or, attrib=value.   More than on attrib is allowed, as a comma-
      // delimited list

      // If same attrib C=2,C=1,D=1
      // e.g. ( C==2 OR C==1 ) AND (D==1)

      NList tlist(0);
      vector<string> ids = tlist.deparseStringList( par::idhelp_lookup_string );
      for ( int i = 0 ; i < ids.size() ; i++)
	{
	  string s = ids[i];
	  if ( s.find("=") == string::npos )
	    error("Lookup query  must be in form ID=value or ATTRIB=value,ATTRIB=value");
	  string fs = s.substr(0, s.find("="));
	  string vs = s.substr(s.find("=")+1);
	  
	  // Does this field exist?
	  map<string,IDField*>::iterator f = fieldMap.find ( fs );
	  if ( f == fieldMap.end() )
	    error("Cannot find field " + fs );
	  
	  IDField * thisField = f->second;

	  // Is this an attribute or an ID?
	  if ( ( ! thisField->attribute ) && ids.size() > 1 )
	    error("Cannot specify multiple ID matches with lookup");
	  
	  IDValue t;
	  t.field = thisField;
	  t.value = vs;	  
	  if (t.field->equiv)
	    t.updateAlias();

	  map<IDField*,set<IDValue> >::iterator i = lookupValues.find( thisField );
	  if ( i == lookupValues.end() )
	    {
	      set<IDValue> ts;
	      ts.insert(t);
	      lookupValues.insert(make_pair(thisField,ts) );
	    }
	  else
	    i->second.insert(t);
	}


      printLOG("Lookup up items matching: ");
      // These will be sorted in field name order, so we
      // can easily figure out OR versus AND conditions
      map<IDField*,set<IDValue> >::iterator i = lookupValues.begin();
      while ( i != lookupValues.end() )
	{
	  set<IDValue>::iterator j = i->second.begin();
	  printLOG( "\n  " + i->first->name + " = " );
	  while ( j != i->second.end() )
	    {
	      printLOG( j->value + " " );
	      ++j;
	    }
	  if ( i->first->attribute ) printLOG(" (attribute)");
	  else printLOG(" (id)");
	  ++i;
	}
      printLOG("\n");
    }
  
  
  printLOG("Writing output to [ " + par::output_file_name + ".id ]\n");
  ofstream OFILE( (par::output_file_name+".id").c_str() , ios::out );

  // Header row
  set<IDField>::iterator f = fields.begin();
  
  while ( f != fields.end() )
    {
      if ( (! par::idhelp_subset ) || subsetFields.find( f->name )!=subsetFields.end() )
	OFILE << f->name << "\t";
      ++f;
    }
  OFILE << "\n";

  
  for ( int g = 0 ; g < idgroup.size(); g++ )
    {
      
      IDGroup * group = &idgroup[g];

      // Has this group already been assigned to a person?
      if ( group->resolved )
	continue;
      
      
      // Make a set of the key groups IDValues
      // If doing a lookup, we might also need to include
      // all fields here

      map<IDField*,IDValue*> keyValues;
      
      for ( int j = 0; j < group->values.size(); j++ )
	{
	  if ( par::idhelp_lookup || 
	       ( ! par::idhelp_subset ) || 
	       subsetFields.find( group->values[j]->field->name )!=subsetFields.end() )
	    keyValues.insert( make_pair( group->values[j]->field, group->values[j] ));
	}


      // If we are filtering, does this person match?

      if ( par::idhelp_lookup )
	{
	  bool match = true;
	  
	  // Compare
	  //  map<IDField*, set<IDValue> > lookupValues;
	  //  with this group
	  
	  map<IDField*, set<IDValue> >::iterator i = lookupValues.begin();

	  while ( i != lookupValues.end() )
	    {

	      map<IDField*,IDValue*>::iterator k = keyValues.find( i->first );
	      if ( k == keyValues.end() ) 
		{
		  match = false;
		}
	      else
		{
		  IDValue * myValue = k->second;
		  
		  set<IDValue>::iterator j = i->second.begin();
		  bool matchField = false;
		  while ( j != i->second.end() )
		    {
		      if ( *myValue == *j )
			matchField = true;
		      ++j;
		    }
		  if ( ! matchField ) 
		    match = false;
		}
	      ++i;
	    }
	  
	  if ( ! match ) 
	    continue;
	}
      
  


      ///////////////////////////////////////////
      // Print row, in same order for all fields

      set<IDField>::iterator f = fields.begin();
      while ( f != fields.end() )
	{
	  
	  if ( (! par::idhelp_subset) || subsetFields.find( f->name )!=subsetFields.end() )	  	  
	    {
	      map<IDField*,IDValue*>::iterator k = keyValues.find( (IDField*)&(*f) ); 
	      if ( k == keyValues.end() )
		OFILE << "." << "\t";
	      else
		OFILE << k->second->value 
		      << " (" 
		      << k->second->jointValue << ")\t";		  
		  
	    }
	  ++f;	    
	}
      OFILE << "\n";
    }  
  OFILE.close();

}
